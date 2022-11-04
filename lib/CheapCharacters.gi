# CheapCharacters.gi        The InduceReduce package         Ulli Kehrle

# these pPrime decomposition functions have been
# implemented by Frank Lübeck
InstallMethod( pPrimepDecompositionMaps,
               [ IsCharacterTable ],

function(tbl)
    return rec ();
end);

InstallMethod(pPrimepDecompositionMap,
              [IsGroup, IsPosInt],

function(G, p)
    local tbl;
    tbl := CharacterTable(G);
    return pPrimepDecompositionMap( tbl, p );
end);


# Let p be a prime.
# An element x of a finite group has a unique decomposition
# x = y z = z y with y being of p'-order and z being a p-element.
# Both, y and z are powers of x, this function returns [k, l] with
# y = x^k and z = x^l. This only depends on the order o of x.
InstallGlobalFunction(pPrimepDecompositionPowers, function(o, p)
  local q, m, g, res;
  if o mod p <> 0 then
    return [1, 0];
  fi;
  q := p;
  m := o/p;
  while m mod p = 0 do
    q := q*p;
    m := m/p;
  od;
  if m = 1 then
    return [0, 1];
  fi;
  g := GcdRepresentation(q, m);
  res := [q*g[1] mod o, m*g[2] mod o];
  if 2*res[1] > o then
    res[1] := res[1]-o;
  fi;
  if 2*res[2] > o then
    res[2] := res[2]-o;
  fi;
  return res;
end);

InstallMethod(pPrimepDecompositionMap,
              [IsCharacterTable, IsPosInt],

function(tbl, p)
    local res, map, es, ords, i;
    map := pPrimepDecompositionMaps(tbl);
    if not IsBound(map.(p)) then
        res := [];
        ords := OrdersClassRepresentatives(tbl);
        for i in [1..Length(ords)] do
            es := pPrimepDecompositionPowers(ords[i], p);
            Add(res, PowerMap(tbl, es[1])[i]);
        od;
        map.(p) := res;
    fi;
    return map.(p);
end);

InstallGlobalFunction(pPrimeCharacter,
function(chi, p)
    local tbl;
    tbl := UnderlyingCharacterTable(chi);
    return Character(tbl, chi{pPrimepDecompositionMap(tbl, p)});
end);

InstallGlobalFunction(pPrimeRestriction,
function(chi, p)
  local s, q, ch, ords, i, tbl;

  tbl := UnderlyingCharacterTable(chi);
  s := Size(tbl);
  q := PPart(s, p);
  ch := [];
  ords := OrdersClassRepresentatives(tbl);
  for i in [1..Length(ords)] do
    if ords[i] mod p = 0 then
      Add(ch, 0);
    else
      Add(ch, q*chi[i]);
    fi;
  od;
  if q*chi = ch then
    ch := chi;
  fi;
  return Character(tbl, ch);
end);

# Code by Frank Lübeck
InstallGlobalFunction(OrderCharacters,
function(GT)
  local ords, sz, fac, e, res, ch, m, val, p, o, k, a;
  if IsOrdinaryTable(GT) then
    ords := OrdersClassRepresentatives(GT);
  elif IsGroup(GT) then
    ords := List(ConjugacyClasses(GT), c-> Order(Representative(c)));
  fi;
  sz := Size(GT);
  fac := Collected(Factors(sz));
  # exponent
  e := Lcm(ords);
  res := [];
  # one class function for each element order
  for o in Set(ords) do
    if o <> 1 then
      ch := [];
      m := e/o;
      val := sz/m;
      for k in ords do
        if m mod k = 0 then
          Add(ch, val);
        else
          Add(ch, 0);
        fi;
      od;
      Add(res, ch);
    fi;
  od;
  return List(res, ch-> Character(GT, ch));
end);


InstallGlobalFunction(IRCharacterToList,
function(GR, chi)
   return Permuted(List(chi), GR!.perm);
end);

InstallGlobalFunction(IRListToCharacter,
function(GR, chi)
   return Character(GR!.C, Permuted(chi, Inverse(GR!.perm)));
end);

InstallGlobalFunction(IRAddIrrChar,
function(GR, char)
  Add(GR!.Ir, IRCharacterToList(char));
end);

InstallGlobalFunction(IRAddChar,
function(GR, char)
  Add(GR!.B, IRCharacterToList(char));
end);

InstallGlobalFunction(IRGenerateMoreCharacters,
function(GR, chars)
    local res, chi, p, i, tbl;

    res := [];
    tbl := GR!.C;

    UniteSet(res, SymmetricParts(tbl, chars, 2));
    UniteSet(res, AntiSymmetricParts(tbl, chars, 2));
    UniteSet(res, SymmetricParts(tbl, chars, 3));
    UniteSet(res, AntiSymmetricParts(tbl, chars, 3));
    UniteSet(res, SymmetricParts(tbl, chars, 4));
    UniteSet(res, AntiSymmetricParts(tbl, chars, 4));


    for chi in chars do
        for p in GR!.Primes do
            AddSet(res, pPrimeCharacter(chi, p));
            AddSet(res, pPrimeRestriction(chi, p));
        od;
    od;


    return res;
end);

InstallGlobalFunction(IRInitStandardChars, function(GR)
    local chars, char, i, m, e, q, n, val, ords, ch, k, o;

    chars := [];

    GR!.Ir := List(LinearCharacters(GR!.G), chi ->  IRCharacterToList(GR, chi));

    # and the regular character;
   # IRAddChar(Order(GR!.G)*IdentityMat(GR!.k)[1]);
    #

    AddSet(chars, TrivialCharacter(GR!.C));

    UniteSet(chars, OrderCharacters(GR!.C));

   if IsPermGroup(GR!.G) then
      # permutation character
      AddSet(chars, NaturalCharacter(GR!.G) );

      # permutation character on pairs
      #if LargestMovedPoint(GR!.G) < 1000 then
      #   AddSet(chars, PermutationCharacter( GR!.G, Tuples(MovedPoints(GR!.G), 2), OnPairs ) );
      #fi;

   elif IsMatrixGroup(GR!.G) then
      #### this is probably not faster than the generic code, even in the presence of a nice
      #### monomorphism

      # if HasNiceMonomorphism( GR!.G ) and IsPermutationGroup(NiceObject( GR!.G )) then
      #    f := NiceMonomorphism( GR!.G );
      #    IRAddChar( List(GR!.reps, x ->
        #         Number( MovedPoints(Image(f)), a -> a = a^(x^f) ) ) );
      # fi;

      # permutation character (OnPoints)
      char := [];
      q := Order( FieldOfMatrixGroup(GR!.G) );
      n := DimensionOfMatrixGroup(GR!.G);
      for x in GR!.classreps do;
         Add(char, q^(n-RankMat(x-IdentityMat(n))));
      od;
      Add(GR!.B, char);

      # permutation character (on 1-dim subspaces)
      char := [];
      for x in GR!.classreps do;
         val := 0;
         for i in Eigenvalues(FieldOfMatrixGroup(GR!.G), x) do
            val := val + q^(n-RankMat(x-i*IdentityMat(n)));
         od;
         Add(char, val);
      od;
      AddSet(chars, char);
   fi;

    UniteSet(chars, IRGenerateMoreCharacters (GR, [TrivialCharacter(GR!.C)]));

    GR!.B := List(chars, chi -> IRCharacterToList(GR, chi));
    IRReduce(GR, [1..GR!.k]);
end);
