#
# InduceReduce.gi        The InduceReduce package         Jonathan Gruber
#
# Implementations
#
# Copyright (C) 2018     Jonathan Gruber
# Copyright (C) 2022     Ulli Kehrle
#
#
IRDefaultOptions := rec(
	IRDoCyclicFirst := true,
    IRDoCyclicLast := false,
    IRLLLOffset := 0,
    IRDelta := 3/4,
    IRIrrAlgorithm := "baum-clausen",
    IRPrimeSelectionStrategy := "exponent",
    IRFusionJobs := 16,
    IRParallelGroups := 8,
);

InstallGlobalFunction(IRGetOpt,
function(optname)
       local val;

       val := ValueOption(optname);
       if val <> fail then
           return val;
       fi;

       return IRDefaultOptions.(optname);
end);

InstallGlobalFunction(IRGetOptions, function()
    local opt, tmp, res;

    res := ShallowCopy(IRDefaultOptions);

    for opt in NamesOfComponents(res) do
        tmp := ValueOption(opt);
        if tmp <> fail then
            res.(opt) := tmp;
        fi;
    od;

    return res;
end);



#############################################################################
##
#F IRInit( <G> )
##
## Initialization of the algorithm, returns a record which contains all
## relevant data concerning the group G
##
InstallGlobalFunction(IRInit,

function(G)
    local GR, PR, i;

    GR := rec();

    GR.G := G;                                             # The group itself
    GR.C := CharacterTable(G);                             # the character table (without irreducible characters)
    GR.n := Size(G);                                       # the order of G
    GR.classes := ShallowCopy(ConjugacyClasses(GR.C));     # the conjugacy classes (mutable)
    GR.k := Size(GR.classes);                              # number of conjugacy classes
    GR.classreps := List(GR.classes,x->Representative(x)); # class representatives
    GR.orders := List(GR.classreps, x->Order(x));          # element orders
        Info(InfoCTUnger, 1, "Induce/Restrict: group with ",
                Length(GR.orders), " conjugacy classes.");
        Info(InfoCTUnger, 2, "Induce/Restrict: orders of class reps: ",
                GR.orders);


    GR.perm := Sortex(GR.orders,function(x,y) return y<x; end);
                                                           # sort by decreasing order of representatives
    GR.classes := Permuted(GR.classes, GR.perm);            # adjust the order of classes and classreps
    GR.classreps := Permuted(GR.classreps, GR.perm);

    GR.powermaps := List(GR.classes, i -> OnTuples(PowerMap(i), GR.perm));
    GR.inverseClasses := List([1..GR.k-1], i -> GR.powermaps[i][GR.orders[i]-1]);
    Add(GR.inverseClasses, GR.k); # unit element

                                                           # permutation that sorts classes by descending order.
    GR.ccsizes := List(GR.classes, x -> Size(x));             # sizes of conjugacy classes
    GR.OrderPrimes := List(GR.orders, x -> Set(Factors(x)));
                                                           # primes dividing the orders of class representatives
    GR.CentralizerPrimes := List([1..GR.k],
        x-> Filtered( Set(Factors(GR.n/GR.ccsizes[x])) , y -> not y in GR.OrderPrimes[x] )  );
                                                           # primes dividing the order of the centralizer,
                                                           # but not the order of the element
    GR.CentralizerPrimeExponents := List([1..GR.k],
        x -> List(GR.CentralizerPrimes[x], p -> PMultiplicity(GR.n/GR.ccsizes[x], p)));
    GR.InduceCyclic := ListWithIdenticalEntries(GR.k, true);
                                                           # the characters of the corresponding cyclic groups still need to be induced
    GR.NumberOfCyclic := GR.k;                             # number of cyclic groups whose characters need to be induced
    GR.NumberOfPrimes := Sum(GR.CentralizerPrimes, x -> Size(x));
                                                           # number of primes for elementary subgroups not used so far
    GR.e := Exponent(G);                                   # exponent of G
    GR.Primes  := PrimeDivisors(GR.e);                    # primes dividing the group order

    GR.Ir := [ Zero([1..GR.k])+1 ];                        # initialize Irr(G) with the trivial character:
    GR.B := [];                                            # B as an empty list
    GR.Gram := [];                                         # Gram empty

    GR.IndexCyc := 0;                                      # position of the cyclic group curretly used in the list of class representatives
    GR.m := 0;                                             # positions from which on characters are not reduced so far

    return Objectify(IRRecordType, GR);
    # return GR;
end);


InstallMethod( ViewObj, "for IRRecord", [ IsIRRecord ],
function(GR)
    local det;

    if GR!.Gram <> [] then
        det := DeterminantMat(GR!.Gram);
    else
        det := 0;
    fi;

    Print("<InduceReduce record with ",
          Size(GR!.Ir), "/", GR!.k, " irreducible characters",
          " and Determinant ", det, " (", FactorsInt(det), "),",
          " cyclics: ", GR!.NumberOfCyclic,
          ", primes: ", GR!.NumberOfPrimes,
          " >");
end);


InstallMethod( PrintObj, "for IRReCord", [ IsIRRecord ],
function(GR)
    local i, j, det;


    if GR!.Gram <> [] then
        det := DeterminantMat(GR!.Gram);
    else
        det := 0;
    fi;

    Print("<InduceReduce record with ",
          Size(GR!.Ir), "/", GR!.k, " irreducible characters",
          " and Determinant ", det, " (", FactorsInt(det), "),"    ,
          " cyclics: ", GR!.NumberOfCyclic,
          ", primes: ", GR!.NumberOfPrimes,
          ".\n"
         );

    Print("\n\nClass Information:");
    for i in [1..GR!.k] do
        if Length(GR!.CentralizerPrimes[i]) >= 1 or GR!.InduceCyclic[i] then
            Print(
                "class ", i,
                ": Order ", GR!.orders[i],
                ", Size ", GR!.ccsizes[i],
                ", Cyclic unused: ", GR!.InduceCyclic[i],
                ", Available Primes: ");
            for j in [1..Length(GR!.CentralizerPrimes[i])] do
                Print(
                    GR!.CentralizerPrimes[i][j], " (",
                    GR!.CentralizerPrimeExponents[i][j], "), "
                    );
            od;
            Print("\n");
        fi;
    od;
    Print(">\n");
end);


InstallGlobalFunction(IRExcludeCovered,
function(GR, Elementary)
        local p;
        # these are (up to conjugacy) contained in Elementary, so they won't
        # yield new information, skip them.
        if IsBound(Elementary.p) then
            IRExcludeElementary(GR,
                                 Elementary.zNr,
                                 Elementary.p,
                                 PPart(GR!.n/GR!.ccsizes[Elementary.zNr],
                                       Elementary.p));
        else
            # remove elementary groups that are completely contained inside
            # the cyclic group Elementary
            for p in GR!.CentralizerPrimes[Elementary.zNr] do
                IRExcludeElementary(GR,
                                    Elementary.zNr,
                                    p,
                                    PPart(GR!.n/GR!.ccsizes[Elementary.zNr], p));
            od;
        fi;
end);


# remove the cyclic groups generated by elements in class classn
# from consideration
InstallGlobalFunction(IRExcludeCyclic,
function(GR, classn)
    local i, powers;

    powers := Set(GR!.powermap[classn]);
    for i in powers do
        if GR!.InduceCyclic[i] then
            Info(InfoCTUnger, 2,
                "Excluding cyclics generated by ",
                i,
                " from consideration");
            GR!.InduceCyclic[i] := false;
            GR!.NumberOfCyclic := GR!.NumberOfCyclic - 1;
        fi;
    od;
end
);

# remove elementary subgroups with cyclic groups inside <classn>
# and p-group of order ppart from consideration
InstallGlobalFunction(IRExcludeElementary,
function(GR, classn, p, ppart)
    local i, pos, powers;
    powers := Set(GR!.powermap[classn]);
    for i in powers do
            if p in GR!.CentralizerPrimes[i] and
                    PPart(GR!.n/ GR!.ccsizes[i], p) = ppart then
                Info(InfoCTUnger, 2,
                    "Excluding elementaries with cyclic part generated by ",
                    i,
                    " and p-group of order ",
                    ppart,
                    " (prime ",
                    p,
                    ") from consideration");

                pos := PositionSorted(GR!.CentralizerPrimes[i], p);
                Remove(GR!.CentralizerPrimes[i], pos);
                Remove(GR!.CentralizerPrimeExponents[i], pos);

                GR!.NumberOfPrimes := GR!.NumberOfPrimes - 1;
            fi;
        od;
    end
);

InstallGlobalFunction( IRClassFusion,
function(GR, Elementary)
    if IsPackageLoaded("IO") then
        return ParListByFork( [1..Elementary.k],
                i -> PositionClass(GR!.G, Elementary.classreps[i],1)^GR!.perm,
                rec(NumberJobs:=16));
    else
        return List( [1..Elementary.k],
                i -> PositionClass(GR!.G, Elementary.classreps[i],1)^GR!.perm);
    fi;
end);

## Compute character table of cyclic group as matrix
InstallGlobalFunction(IRVandermonde , function(GR, Elementary)
    local i, j, M;
        M:=NullMat(GR!.orders[Elementary.zNr],GR!.orders[Elementary.zNr]);
        for i in [0..GR!.orders[Elementary.zNr]-1] do
            for j in [0..GR!.orders[Elementary.zNr]-1] do
                M[i+1][j+1]:=E(GR!.orders[Elementary.zNr])^(i*j);
            od;
        od;
        return M;
    end);

# conjugacy class representatives of cyclic group corresponding to
# the ordering of columns in Vandermonde
InstallGlobalFunction(IRClassesCyclic, function(GR, Elementary)
    local res, h;
        res:=[Identity(GR!.G)];
        h:=Elementary.z;
        while not h=Identity(GR!.G) do
            Add(res,h);
            h:=h*Elementary.z;
        od;
        return res;
    end
);


InstallGlobalFunction(IRNewElementary ,

function(GR, classn, p)
    local Elementary, pos;

    Info(InfoCTUnger, 2,
            "Considering new elementary group with prime ", p,
            " and cyclic generator class ", classn);

    Elementary := rec ();

    Elementary.isCyclic := false;
    Elementary.zNr := classn;
    Elementary.z := GR!.classreps[classn];
    Elementary.p := p;

    GR!.NumberOfPrimes := GR!.NumberOfPrimes - 1;
    pos := PositionSorted(GR!.CentralizerPrimes[classn] , Elementary.p);

    Remove(GR!.CentralizerPrimes[classn], pos);
    Remove(GR!.CentralizerPrimeExponents[classn], pos);

    return Elementary;
end);

InstallGlobalFunction(IRNewCyclic,

function(GR, classn)
    local Elementary;
    Elementary := rec ();

    Elementary.isCyclic := true;
    Elementary.zNr := classn;
    Elementary.z := GR!.classreps[classn];
    GR!.InduceCyclic[classn] := false;
    GR!.NumberOfCyclic := GR!.NumberOfCyclic - 1;

    return Elementary;
end);


# IR.FindGroup
# constructs the record for the next group to induce from
InstallGlobalFunction(    IRFindGroup ,

function(GR)
    local p, k, Opt, primes, i, det, classes;

    Opt := IRGetOptions();

    if GR!.NumberOfCyclic = 0 and GR!.NumberOfPrimes = 0 then
        return fail;
    fi;

    if GR!.NumberOfPrimes <> 0 then
        if Size(GR!.Ir) + Size(GR!.B) = GR!.k then
            det := DeterminantMat(GR!.Gram);

            # primes where we still need something.
            primes := Filtered(GR!.Primes, p -> det mod p = 0 );

            classes := Filtered([1..GR!.k], i -> Intersection(GR!.CentralizerPrimes[i], primes) <> [] );

            if classes <> [] then
                i := PseudoRandom(classes);
                p := PseudoRandom(Intersection(GR!.CentralizerPrimes[i], primes));

                return IRNewElementary(GR, i, p);
            fi;
        fi;
    fi;

    while true do
        # run over the conjugacy classes
        GR!.IndexCyc := GR!.IndexCyc + 1;
        if GR!.IndexCyc > GR!.k then
            GR!.IndexCyc := 1;
        fi;

        # the class to use for the cyclic group
        k := GR!.IndexCyc;

        if not IsEmpty(GR!.CentralizerPrimes[k]) then
            if Opt.IRPrimeSelectionStrategy = "random" then
                p := PseudoRandom(GR!.CentralizerPrimes[k]);
                return IRNewElementary(GR, k, p);
            elif Opt.IRPrimeSelectionStrategy = "exponent" then
                p := GR!.CentralizerPrimes[k][PositionMaximum(GR!.CentralizerPrimeExponents[k])];
                return IRNewElementary(GR, k, p);
            else
                Error("Unknown strategy");
            fi;
        elif GR!.InduceCyclic[k] and
            (not Opt.IRDoCyclicLast or GR!.NumberOfPrimes<=0) then
            # if there are no primes left for this class, induce from cyclic group
            # if Opt.DoCyclicLast=true: only induce from cyclic groups
            # when all primes have been used
            return IRNewCyclic(GR, k);
        fi;
    od;
end);



# returns record with components:
# - characters
# - fusedClasses
InstallGlobalFunction(  IRHandleElementary ,

function(GR, Elementary)
    local i,j,i1,j1,p,temp,powermap,mat;

    if Elementary.isCyclic then
        Elementary.n:=GR!.orders[Elementary.zNr]; # order of the elementary group
        Elementary.k:=Elementary.n; # number of conjugacy classes
        Elementary.ccsizes:=ListWithIdenticalEntries(Elementary.k,1); # class sizes
        powermap:=OnTuples(Concatenation([1],PowerMap(GR!.classes[Elementary.zNr])), GR!.perm);
        Elementary.classfusion:=powermap;
            # class fusion equals power map for cyclic group
        Elementary.classreps:=IRClassesCyclic(GR, Elementary); # class representatives
        Elementary.XE:=IRVandermonde(GR, Elementary); # character table
    else
        Elementary.P := SylowSubgroup(Centralizer(GR!.G, GR!.classreps[Elementary.zNr]), Elementary.p);
        Elementary.n:=GR!.orders[Elementary.zNr]*Size(Elementary.P); # order
        Elementary.ctblP:=CharacterTable(Elementary.P);
            # character table of p-group
        Elementary.classrepsP:=List(ConjugacyClasses(Elementary.ctblP),
            x->Representative(x));
            #classes of p-group in corresponding order
        Elementary.ccsizesP:=List(ConjugacyClasses(Elementary.ctblP),x->Size(x));                 # class sizes p-group
        Elementary.kP:=Size(Elementary.classrepsP); # number of classes of p-group
        Elementary.classrepsZ:=IRClassesCyclic(GR, Elementary);
            # class representatives cyclic group
        Elementary.kZ:=GR!.orders[Elementary.zNr]; # number of classes cyclic group
        Elementary.k:=Elementary.kP*Elementary.kZ;
                                # number of classes elementary group
        Elementary.classreps := []; # compute the class representatives
        Elementary.ccsizes := []; # and the class sizes
        for i in [1..Elementary.kP] do
            for j in [1..Elementary.kZ] do
                Add(Elementary.classreps ,
                    Elementary.classrepsP[i]*Elementary.classrepsZ[j]);
                Add(Elementary.ccsizes,Elementary.ccsizesP[i]);
            od;
        od;
        Elementary.classfusion := IRClassFusion(GR, Elementary);
        # f := IsomorphismPcGroup(Elementary.P)
        Elementary.XP := IrrBaumClausen(Elementary.P);
        Elementary.XZ := IRVandermonde(GR, Elementary); # character table of the cycic group
        Elementary.XE := []; # compute character table of the elementary group
        for i in [1..Elementary.kP] do
            for j in [1..Elementary.kZ] do
                temp:=[];
                for i1 in [1..Elementary.kP] do
                    for j1 in [1..Elementary.kZ] do
                        Add(temp,Elementary.XP[i][i1]*Elementary.XZ[j][j1]);
                    od;
                od;
                Add(Elementary.XE,temp);
            od;
        od;
    fi;

    mat:=NullMat(Elementary.k,GR!.k);
    for i in [1..Elementary.k] do
        for j in [1..Elementary.k] do
            mat[i][Elementary.classfusion[j]]:=mat[i][Elementary.classfusion[j]]+
                ( (GR!.n/GR!.ccsizes[Elementary.classfusion[j]]) /
                (Elementary.n/Elementary.ccsizes[j])*Elementary.XE[i][j] );
        od;
    od;
    mat := Set(mat); # remove duplicates

    return rec ( characters := mat, fusedClasses := Set(Elementary.classfusion) );
end);


#############################################################################
##
#F IRDoCyclics( <GR> )
##
## induce characters from all cyclic subgroups
##
InstallGlobalFunction( IRDoCyclics,

function(GR)
    local ords, inds;
    ords := ShallowCopy(OrdersClassRepresentatives(GR!.C));
    inds := [1..NrConjugacyClasses(GR!.C)];
    SortParallel(ords, inds, function(x,y) return x>y; end);
    Append(GR!.B, List(InducedCyclic(GR!.C, inds, "all"), ch -> Permuted(ch, GR!.perm)));

    GR!.InduceCyclic:=ListWithIdenticalEntries(GR!.k,false);
    GR!.NumberOfCyclic:=0;

    IRReduce(GR, [1..GR!.k]);

    return;
end);



#############################################################################
##
#F IR.tinI( <GR> )
##
## adjusts the signs of the characters and permutes them to the right
## ordering of conjugacy classes
##
InstallGlobalFunction(IRtinI, function(GR)
    local irr,i,perm, T;

        for i in [1..GR!.k] do # adjust the signs
            if GR!.Ir[i][GR!.k]<0 then
                GR!.Ir[i]:=-GR!.Ir[i];
            fi;
        od;
        irr:=[];
        for i in [1..GR!.k] do
            # permute irreducible characters back to the order of classes in GR.C
            GR!.Ir[i]:=Permuted(GR!.Ir[i],Inverse(GR!.perm));
            irr[i]:=Character(GR!.C,GR!.Ir[i]);
        od;


        # convert the result into a chacacter table
        T:=rec();
        T.UnderlyingCharacteristic:=0;
        T.NrConjugacyClasses:=GR!.k;
        T.ConjugacyClasses:=ConjugacyClasses(GR!.C);
        T.Size:=GR!.n;
        T.OrdersClassRepresentatives:=Permuted(GR!.orders,Inverse(GR!.perm));
        T.SizesCentralizers:=Permuted(List(GR!.ccsizes,x->GR!.n/x),Inverse(GR!.perm));
        T.Irr:=irr;
        T.UnderlyingGroup:=UnderlyingGroup(GR!.C);
        T.IdentificationOfConjugacyClasses:=IdentificationOfConjugacyClasses(GR!.C);
        T.InfoText:="Computed using Unger's induce-reduce algorithm";
        ConvertToCharacterTableNC(T);
        return T;
end );

InstallGlobalFunction (IRIsComplete,
function(GR)
  local det;

    # if GR!.Gram <> [] then
    #     det := DeterminantMat(*GR!.Gram);
    # else
    #     det := 0;
    # fi;

    Info(InfoCTUnger, 2, "Progress: ", GR!.NumberOfCyclic, " / ", GR!.NumberOfPrimes, "\n");

    # Info(InfoCTUnger, 2, "Checking Progress of ", GR,
    #      ". Ir: ", Size(GR!.Ir),
    #      ", B: ", Size(GR!.B),
    #      ", m: ", GR!.m,
    #      ", det: ", det,
    #      "\n"
    #     );

    # we're done if we have all the irreducible characters or no more groups
    return Size(GR!.Ir) >= GR!.k;
    #or (GR!.NumberOfPrimes = 0 and GR!.NumberOfCyclic = 0);

    # alternative approach: we're done if we have found the whole character ring
    # return Size(GR!.Ir) + Size(GR!.B) >= GR!.k and det in [0,1];
    end);


#############################################################################
##
#F InduceReduce(GR)
##
## computes the irreducible characters of the group GR!.G (upto sign) and puts
## them in the component GR!.Ir of GR, returns GR!.Ir
##
InstallGlobalFunction( InduceReduce,
function(GR)
    local H, ccsizesH, temp, it, Elementary, tmp, p, i, Opt;

    Opt := IRGetOptions();

    if Opt.IRDoCyclicFirst = true then
        Info(InfoCTUnger, 1, "Induce: from cyclic subgroups");
        IRDoCyclics(GR);
    fi;

    if Size(GR!.Ir) >= GR!.k then
        return GR!.Ir;
    fi;

    while not IRIsComplete(GR) do
        # Opt.LLLOffset postpones the first LLL lattice reduction
        Opt.IRLLLOffset := Opt.IRLLLOffset - 1;

        Elementary := IRFindGroup(GR); # find elementary subgroup

        IRExcludeCovered(GR, Elementary);

        tmp := IRHandleElementary(GR, Elementary);

        # add our new characters
        Append(GR!.B, tmp.characters);

        Info(InfoCTUnger, 1, "Induce/Restrict: Trying [|Z|, |P|, k(E)] = ",
            [ GR!.orders[GR!.IndexCyc],
                Elementary.n/GR!.orders[GR!.IndexCyc], Elementary.k ]);

        IRReduce(GR, tmp.fusedClasses : Opt); # reduce GR!.B by GR!.Ir and do lattice reduction
    od;

    if Size(GR!.Ir)>=GR!.k then
        return GR!.Ir;
    else
        # if there are still characters missing, try one last LLL reduction with Opt.Delta = 1
        Opt.IRDelta:=1;
        Opt.IRLLLOffset := 0;
        IRReduce(GR, tmp.fusedClasses : Opt);
        if Size(GR!.Ir)>=GR!.k then
            return GR!.Ir;
        fi;

        OnBreak := function() end;
        OnBreakMessage := function() end;
        Error("**********************************************************",
              "\n\n",
              "Run out of subgroups and still did not find all irreducible ",
              "characters!\n",
              "They are guaranteed to be integer linear combinations of\n",
              "GR!.Ir and GR!.B, but lattice reduction did not find them.\n",
              "You now have the opportunity to alter GR!.Ir in the REPL.\n",
              "Type \"return;\" when you're done.\n\n");

        return GR!.Ir;
    fi;

end );


#############################################################################
##
#F CharacterTableUnger( <G> )
## 
## Computes the character table of a finite group using Unger's algorithm
##
InstallGlobalFunction( CharacterTableUnger,
function(G, Options...)
local GR, T;

    GR := IRInit(G);
    IRInitStandardChars(GR);

    if IsPackageLoaded("IO") then
        InduceReduceParallel(GR);
    else
        InduceReduce(GR); # do the induce-reduce algorithm
    fi;

	return IRtinI(GR);
end );
