################################################################################
##
#F  BCLInit( <G> )  . . . . . . . . . inits record for layered Baum-Clausen algo
##
InstallGlobalFunction( BCLInit,

function( G )
    local e,             # occurring roots of unity are `e'-th roots
          pcgs,          # Pcgs of `G'
          abel,          # position of abelian normal comp. subgroup
          ExtLinRep,     # local function
          linear,        # list of linear representations
          p,             # size of current composition factor
          pexp,          # exponent vector of `pcgs[i]^p'
          root,          # value of an extension
          roots,         # list of $p$-th roots (relative to `e')
          mulmoma,       # product of two monomial matrices
          poweval,       # representing matrix for power of generator
          i, d, j, k, l, # loop variables
          M,             #
          ssr,           #

          BCR;           # the record holding all the information

    # create the record
    BCR := rec();

    BCR.G := G;
    BCR.complete := false;

    if not IsSolvableGroup( G ) then
      Error( "<G> must be solvable" );
    fi;


    # Step 0:
    # Treat the case of the trivial group,
    # and initialize some variables.

    BCR.pcgs:= SpecialPcgs( G );
#T because I need a ``prime orders pcgs''
    BCR.lg:= Length( BCR.pcgs );

    if BCR.lg = 0 then
        BCR.kernel := G;
        BCR.exponent := 1;
        BCR.nonlin := [];
        BCR.lin := [ [] ];

        BCR.complete := true;
        return BCR;
    fi;

    BCR.cs:= PcSeries( BCR.pcgs );

    if HasExponent( G ) then
      BCR.exponent:= Exponent( G );
    else
      BCR.exponent:= Size(G);
#T better adjust on the fly
    fi;


    # Step 1:
    # If necessary then adjust the composition series of $G$
    # and get the largest factor group of $G$ that has an abelian normal
    # subgroup such that the factor group modulo this subgroup is
    # supersolvable.

    abel:= 1;
    while IsNormal( G, BCR.cs[ abel ] ) and not IsAbelian( BCR.cs[ abel ] ) do
      abel:= abel + 1;
    od;

    BCR.abel := abel;

    # If `cs[ abel ]' is abelian then we compute its representations first,
    # and then loop over the initial part of the composition series;
    # note that the factor group is supersolvable.
    # If `cs[ abel ]' is not abelian then we try to switch to a better
    # composition series, namely one through the derived subgroup of the
    # supersolvable residuum.

    if not IsNormal( G, BCR.cs[ abel ] ) then

      # We have reached a non-normal nonabelian composition subgroup
      # so we have to adjust the composition series.

      Info( InfoGroup, 2,
            "BaumClausenInfo: switching to a suitable comp. ser." );

      ssr:= SupersolvableResiduumDefault( G );
      BCR.hom:= NaturalHomomorphismByNormalSubgroupNC( G,
                DerivedSubgroup( ssr.ssr ) );

      # `SupersolvableResiduumDefault' contains a component `ds',
      # a list of subgroups such that any composition series through
      # `ds' from `G' down to the residuum is a chief series.
      pcgs:= [];
      for i in [ 2 .. Length( ssr.ds ) ] do
        j:= NaturalHomomorphismByNormalSubgroupNC( ssr.ds[ i-1 ], ssr.ds[i] );
        Append( pcgs, List( SpecialPcgs( ImagesSource( j ) ),
                            x -> PreImagesRepresentative( j, x ) ) );
      od;
      Append( pcgs, SpecialPcgs( ssr.ds[ Length( ssr.ds ) ]) );
      G:= ImagesSource( BCR.hom );

      BCR.pcgs:= List( pcgs, x -> ImagesRepresentative( BCR.hom, x ) );
      BCR.pcgs:= Filtered( BCR.pcgs, x -> Order( x ) <> 1 );
      BCR.pcgs:= PcgsByPcSequence( ElementsFamily( FamilyObj( G ) ), BCR.pcgs );
      BCR.cs:= PcSeries( BCR.pcgs );
      BCR.lg:= Length( BCR.pcgs );

      # The image of `ssr' under `hom' is abelian,
      # compute its position in the composition series.
      BCR.abel:= Position( BCR.cs, ImagesSet( BCR.hom, ssr.ssr ) );

      # If `G' is supersolvable then `abel = lg+1',
      # but the last *nontrivial* member of the chain is normal and abelian,
      # so we choose this group.
      # (Otherwise we would have the technical problem in step 4 that the
      # matrix `M' would be empty.)
      if BCR.lg < BCR.abel then
        BCR.abel:= BCR.lg;
      fi;

    fi;

    # Step 2:
    # Compute the representations of `cs[ abel ]',
    # each a list of images of $g_{abel}, \ldots, g_{lg}$.

    # The local function `ExtLinRep' computes the extensions of the
    # linear $G_{i+1}$-representations $F$ in the list `linear' to $G_i$.
    # The only condition that must be satisfied is that
    # $F(g_i)^p = F(g_i^p)$.
    # (Roughly speaking, we just compute $p$-th roots.)

    ExtLinRep:= function( i, linear, pexp, roots )

      local nextlinear, rep, j, shift;

      nextlinear:= [];
      if IsZero( pexp ) then

        # $g_i^p$ is the identity
        for rep in linear do
          for j in roots do
            Add( nextlinear, Concatenation( [ j ], rep ) );
          od;
        od;

      else

        pexp:= pexp{ [ i+1 .. BCR.lg ] };
#T cut this outside the function!
        for rep in linear do

          # Compute the value of `rep' on $g_i^p$.
          shift:= pexp * rep;

          if shift mod p <> 0 then
            # We must enlarge the exponent.
            Error("wrong exponent");
#T if not integral then enlarge the exponent!
#T (is this possible here at all?)
          fi;
          shift:= shift / p;
          for j in roots do
            Add( nextlinear, Concatenation( [ (j+shift) mod BCR.exponent ], rep ) );
          od;

        od;

      fi;

      return nextlinear;
    end;


    BCR.indices:= RelativeOrders( BCR.pcgs );
#T here set the exponent `e' to `indices[ lg ]' !
    Info( InfoGroup, 2,
          "BaumClausenInfo: There are ", BCR.lg, " steps" );

    BCR.linear:= List( [ 0 .. BCR.indices[BCR.lg]-1 ] * ( BCR.exponent / BCR.indices[BCR.lg] ),
                   x -> [ x ] );


    for i in [ BCR.lg-1, BCR.lg-2 .. BCR.abel ] do

      Info( InfoGroup, 2,
            "BaumClausenInfo: Compute repres. of step ", i,
            " (central case)" );

      p:= BCR.indices[i];

      # `pexp' describes $g_i^p$.
      pexp:= ExponentsOfRelativePower( BCR.pcgs,i);
# { ? } ??

      root:= BCR.exponent/p;
#T enlarge the exponent if necessary!
      roots:= [ 0, root .. (p-1)*root ];
      BCR.linear:= ExtLinRep( i, BCR.linear, pexp, roots );

    od;

    # We are done if $G$ is abelian.
    if BCR.abel = 1 then
        BCR.kernel := TrivialSubgroup( G );
        BCR.exponent := e;
        BCR.nonlin := [];
        BCR.lin := linear;

        BCR.complete := true;
        return BCR;
    fi;


    # Step 3:
    # Define some local functions.
    # (We did not need them for abelian groups.)

    # `mulmoma' returns the product of two monomial matrices.
    mulmoma:= function( a, b )
      local prod, i;
      prod:= rec( perm := b.perm{ a.perm },
                  diag := [] );
      for i in [ 1 .. Length( a.perm ) ] do
        prod.diag[ b.perm[i] ]:= ( b.diag[ b.perm[i] ] + a.diag[i] ) mod e;
      od;
      return prod;
    end;

    # `poweval' evaluates the representation `rep' on the $p$-th power of
    # the conjugating element.
    # This $p$-th power is described by `poli'.
    poweval:= function( rep, poli )
      local pow, i;
      if IsEmpty( poli ) then
        return rec( perm:= [ 1 .. Length( rep[1].perm ) ],
                    diag:= [ 1 .. Length( rep[1].perm ) ] * 0 );
      fi;
      pow:= rep[ poli[1] ];
      for i in [ 2 .. Length( poli ) ] do
        pow:= mulmoma( pow, rep[ poli[i] ] );
      od;
      return pow;
    end;


    # Step 4:
    # Compute the actions of $g_j$, $j < abel$, on the representations
    # of $G_{abel}$.
    # Let $g_i^{g_j} = \prod_{k=1}^n g_k^{\alpha_{ik}^j}$,
    # and set $A_j = [ \alpha_{ik}^j} ]_{i,k}$.
    # Then the representation that maps $g_i$ to the root $\zeta_e^{c_i}$
    # is mapped to the representation that has images exponents
    # $A_j * (c_1, \ldots, c_n)$ under $g_j$.

    Info( InfoGroup, 2,
          "BaumClausenInfo: Initialize actions on abelian normal subgroup" );

    BCR.pilinear:= [];
    for j in [ 1 .. BCR.abel-1 ] do

      # Compute the matrix $A_j$.
      M:= List( [ BCR.abel .. BCR.lg ],
                i -> ExponentsOfPcElement( BCR.pcgs, BCR.pcgs[i]^BCR.pcgs[j],
                                           [ BCR.abel .. BCR.lg ] ) );

      # Compute the permutation corresponding to the action of $g_j$.
      BCR.pilinear[j]:= List( BCR.linear,
                          rep -> Position( BCR.linear,
                                           List( M * rep, x -> x mod BCR.exponent ) ) );

    od;


    BCR.nonlin   := [];
    BCR.pinonlin := [];
    BCR.Xlist    := [];

    BCR.position := BCR.abel-1;

    return BCR;
end);


################################################################################
##
#F  BCLStep( <BCR> )  . . . . . . . . . . . . performs step of Baum-Clausen algo
##
# InstallGlobalFunction( BCLStep,
function(BCR)
    local e,             # occurring roots of unity are `e'-th roots
          pcgs,          # Pcgs of `G'
          lg,            # length of `pcgs'
          linear,        # list of linear representations
          i,             # current position in the iteration: $G_i$
          p,             # size of current composition factor
          pexp,          # exponent vector of `pcgs[i]^p'
          root,          # value of an extension
          roots,         # list of $p$-th roots (relative to `e')
          mulmoma,       # product of two monomial matrices
          poweval,       # representing matrix for power of generator
          pilinear,      # action of $g_1, \ldots, g_i$ on `linear'
          d, j, k, l,    # loop variables
          u, v, w,       # loop variables
          pos,           # position in a list
          nonlin,        # list of nonlinear representations
          pinonlin,      # action of $g_1, \ldots, g_i$ on `nonlin'
          Xlist,         # conjugating matrices:
                         # for $X = `Xlist[j][k]'$, we have
                         # $X \cdot {`nonlin[k]'}^{g_j} \cdot X^{-1} =
                         #    `nonlin[ pinonlin[j][k] ]'$
          min,           #
          minval,        #
          nextlinear,    # extensions of `linear'
          nextnonlin1,   # nonlinear repr. arising from `linear'
          nextnonlin2,   # nonlinear repr. arising from `nonlin'
          pinextlinear,  # action of $g_1, \ldots, g_i$ on `nextlinear'
          pinextnonlin1, # action of $g_1, \ldots, g_i$ on `nextnonlin1'
          pinextnonlin2, # action of $g_1, \ldots, g_i$ on `nextnonlin2'
          nextXlist1,    # conjugating matrices for `nextnonlin1'
          nextXlist2,    # conjugating matrices for `nextnonlin2'
          cexp,          # exponent vector of `pcgs[i]^pcgs[j]'
          poli,          # list that encodes `pexp'
          rep,           # one representation
          D, C,          #
          value,         #
          image,         #
          used,          # boolean list
          Dpos1,         # positions of extension resp. induced repres.
                         # that arise from linear representations
          Dpos2,         # positions of extension resp. induced repres.
                         # that arise from nonlinear representations
          dim,           # dimension of the current representation
          invX,          # inverse of `X'
          D_gi,          #
          orb,           #
          Forb,          #
          sigma, pi,     # permutations needed in the fusion case
          constants,     # vector $(c_0, c_1, \ldots, c_{p-1})$
          kernel;        # kernel of `hom'


    # Step 3:
    # Define some local functions.
    # (We did not need them for abelian groups.)

    # `mulmoma' returns the product of two monomial matrices.
    mulmoma:= function( a, b )
      local prod, i;
      prod:= rec( perm := b.perm{ a.perm },
                  diag := [] );
      for i in [ 1 .. Length( a.perm ) ] do
        prod.diag[ b.perm[i] ]:= ( b.diag[ b.perm[i] ] + a.diag[i] ) mod e;
      od;
      return prod;
    end;

    # `poweval' evaluates the representation `rep' on the $p$-th power of
    # the conjugating element.
    # This $p$-th power is described by `poli'.
    poweval:= function( rep, poli )
      local pow, i;
      if IsEmpty( poli ) then
        return rec( perm:= [ 1 .. Length( rep[1].perm ) ],
                    diag:= [ 1 .. Length( rep[1].perm ) ] * 0 );
      fi;
      pow:= rep[ poli[1] ];
      for i in [ 2 .. Length( poli ) ] do
        pow:= mulmoma( pow, rep[ poli[i] ] );
      od;
      return pow;
    end;




    # Step 5:
    # Run up the composition series from `abel' to `1',
    # and compute extensions resp. induced representations.
    # For each index, we have to update `linear', `pilinear',
    # `nonlin', `pinonlin', and `Xlist'.
    #
    #

    # already finished
    if BCR.complete = true or BCR.position <= 0 then
        return BCR;
    fi;

    linear   := BCR.linear;
    pilinear := BCR.pilinear;
    nonlin   := BCR.nonlin;
    pinonlin := BCR.pinonlin;
    Xlist    := BCR.Xlist;

    i := BCR.position;

    e := BCR.exponent;

    p := BCR.indices[i];
    pcgs := BCR.pcgs;
    lg := BCR.lg;

      # `poli' describes $g_i^p$.
      #was pexp:= ExponentsOfPcElement( pcgs, pcgs[i]^p );
      pexp:= ExponentsOfRelativePower( pcgs, i );
      poli:= Concatenation( List( [ i+1 .. lg ],
                                  x -> List( [ 1 .. pexp[x] ],
                                             y -> x-i ) ) );

      # `p'-th roots of unity
      roots:= [ 0 .. p-1 ] * ( e/p );

      Info( InfoGroup, 2,
            "BaumClausenInfo: Compute repres. of step ", i );

      # Step A:
      # Compute representations of $G_i$ arising from *linear*
      # representations of $G_{i+1}$.

      used        := BlistList( [ 1 .. Length( linear ) ], [] );
      nextlinear  := [];
      nextnonlin1 := [];
      d           := 1;

      pexp:= pexp{ [ i+1 .. lg ] };

      # At position `d', store the position of either the first extension
      # of `linear[d]' in `nextlinear' or the position of the induced
      # representation of `linear[d]' in `nextnonlin1'.
      Dpos1:= [];

      while d <> fail do

        rep:= linear[d];
        used[d]:= true;

        # `root' is the value of `rep' on $g_i^p$.
        root:= ( pexp * rep ) mod e;

        if pilinear[i][d] = d then

          # `linear[d]' extends to $G_i$.
          Dpos1[d]:= Length( nextlinear ) + 1;

          # Take a `p'-th root.
          root:= root / p;
#T enlarge the exponent if necessary!

          for j in roots do
            Add( nextlinear, Concatenation( [ root+j ], rep ) );
          od;

        else

          # We must fuse the representations in the orbit of `d'
          # under `pilinear[i]';
          # so we construct the induced representation `D'.

          Dpos1[d]:= Length( nextnonlin1 ) + 1;

          D:= List( rep, x -> rec( perm := [ 1 .. p ],
                                   diag := [ x ]
                                  ) );
          pos:= d;
          for j in [ 2 .. p ] do

            pos:= pilinear[i][ pos ];
            for k in [ 1 .. Length( rep ) ] do
              D[k].diag[j]:= linear[ pos ][k];
            od;
            used[ pos ]:= true;
            Dpos1[ pos ]:= Length( nextnonlin1 ) + 1;

          od;

          Add( nextnonlin1,
               Concatenation( [ rec( perm := Concatenation( [p], [1..p-1]),
                                     diag := Concatenation( [ 1 .. p-1 ] * 0,
                                                            [ root ] ) ) ],
                              D ) );
          Assert( 2, BaumClausenInfoDebug.testrep( pcgs{ [ i .. lg ] },
                              nextnonlin1[ Length( nextnonlin1 ) ], e ),
                  Concatenation( "BaumClausenInfo: failed assertion in ",
                      "inducing linear representations ",
                      "(i = ", String( i ), ")\n" ) );

        fi;

        d:= Position( used, false, d );

      od;


      # Step B:
      # Now compute representations of $G_i$ arising from *nonlinear*
      # representations of $G_{i+1}$ (if there are some).

      used:= BlistList( [ 1 .. Length( nonlin ) ], [] );
      nextnonlin2:= [];
      if Length( nonlin ) = 0 then
        d:= fail;
      else
        d:= 1;
      fi;

      # At position `d', store the position of the first extension resp.
      # of the induced representation of `nonlin[d]'in `nextnonlin2'.
      Dpos2:= [];

      while d <> fail do

        used[d]:= true;
        rep:= nonlin[d];

        if pinonlin[i][d] = d then

          # The representation $F = `rep'$ has `p' different extensions.
          # For `X = Xlist[i][d]', we have $`rep ^ X' = `rep'^{g_i}$,
          # i.e., $X^{-1} F X = F^{g_i}$.
          # Representing matrix $F(g_i)$ is $c X$ with $c^p X^p = F(g_i^p)$,
          # so $c^p X^p.diag[k] = F(g_i^p).diag[k]$ for all $k$ ;
          # for determination of $c$ we look at `k = X^p.perm[1]'.

          X:= Xlist[i][d];
          image:= X.perm[1];
          value:= X.diag[ image ];
          for j in [ 2 .. p ] do

            image:= X.perm[ image ];
            value:= X.diag[ image ] + value;
            # now `image = X^j.perm[1]', `value = X^j.diag[ image ]'

          od;

          # Subtract this from $F(g_i^p).diag[k]$;
          # note that `image' is the image of 1 under `X^p', so also
          # under $F(g_i^p)$.
          value:= - value;
          image:= 1;
          for j in poli do
            image:= rep[j].perm[ image ];
            value:= rep[j].diag[ image ] + value;
          od;

          value:= ( value / p ) mod e;
#T enlarge the exponent if necessary!

          Dpos2[d]:= Length( nextnonlin2 ) + 1;

          # Compute the `p' extensions.
          for k in roots do
            Add( nextnonlin2, Concatenation(
                    [ rec( perm := X.perm,
                      diag := List( X.diag,
                             x -> ( x  + k + value ) mod e ) ) ], rep ) );
            Assert( 2, BaumClausenInfoDebug.testrep( pcgs{ [ i .. lg ] },
                                nextnonlin2[ Length( nextnonlin2 ) ], e ),
                    Concatenation( "BaumClausenInfo: failed assertion in ",
                        "extending nonlinear representations ",
                        "(i = ", String( i ), ")\n" ) );
          od;

        else

          # `$F$ = nonlin[d]' fuses with `p-1' partners given by the orbit
          # of `d' under `pinonlin[i]'.
          # The new irreducible representation of $G_i$ will be
          # $X Ind( F ) X^{-1}$ with $X$ the block diagonal matrix
          # consisting of blocks $X_{i,F}^{(k)}$ defined by
          # $X_{i,F}^{(0)} = Id$,
          # and $X_{i,F}^{(k)} = X_{i,\pi_i^{k-1} F} X_{i,F}^{(k-1)}$
          # for $k > 0$.

          # The matrix for $g_i$ in the induced representation $Ind( F )$ is
          # of the form
          #       | 0   F(g_i^p) |
          #       | I      0     |
          # Thus $X Ind(F) X^{-1} ( g_i )$ is the block diagonal matrix
          # consisting of the blocks
          # $X_{i,F}, X_{i,\pi_i F}, \ldots, X_{i,\pi_i^{p-2} F}$, and
          # $F(g_i^p) \cdot ( X_{i,F}^{(p-1)} )^{-1}$.

          dim:= Length( rep[1].diag );
          Dpos2[d]:= Length( nextnonlin2 ) + 1;

          # We make a copy of `rep' because we want to change it.
          D:= List( rep, record -> rec( perm := ShallowCopy( record.perm ),
                                        diag := ShallowCopy( record.diag )
                                       ) );

          # matrices for $g_j, i\< j \leq n$
          pos:= d;
          for j in [ 1 .. p-1 ] * dim do
            pos:= pinonlin[i][ pos ];
            for k in [ 1 .. Length( rep ) ] do
              Append( D[k].diag, nonlin[ pos ][k].diag );
              Append( D[k].perm, nonlin[ pos ][k].perm + j );
            od;

            used[ pos ]:= true;
            Dpos2[ pos ]:= Length( nextnonlin2 ) + 1;

          od;

          # The matrix of $g_i$ is a block-cycle with blocks
          # $X_{i,\pi_i^k(F)}$ for $0 \leq k \leq p-2$,
          # and $F(g_i^p) \cdot (X_{i,F}^{(p-1)})^{-1}$.

          X:= Xlist[i][d];      # $X_{i,F}$
          pos:= d;
          for j in [ 1 .. p-2 ] do
            pos:= pinonlin[i][ pos ];
            X:= mulmoma( Xlist[i][ pos ], X );
          od;

          # `invX' is the inverse of `X'.
          invX:= rec( perm := [], diag := [] );
          for j in [ 1 .. Length( X.diag ) ] do
            invX.perm[ X.perm[j] ]:= j;
            invX.diag[j]:= e - X.diag[ X.perm[j] ];
          od;
#T improve this using the {} operator!

          X:= mulmoma( poweval( rep, poli ), invX );
          D_gi:= rec( perm:= List( X.perm, x -> x  + ( p-1 ) * dim ),
                      diag:= [] );

          pos:= d;
          for j in [ 0 .. p-2 ] * dim do

            # $X_{i,\pi_i^j F}$
            Append( D_gi.diag, Xlist[i][ pos ].diag);
            Append( D_gi.perm, Xlist[i][ pos ].perm + j);
            pos:= pinonlin[i][ pos ];

          od;

          Append( D_gi.diag, X.diag );

          Add( nextnonlin2, Concatenation( [ D_gi ], D ) );
          Assert( 2, BaumClausenInfoDebug.testrep( pcgs{ [ i .. lg ] },
                              nextnonlin2[ Length( nextnonlin2 ) ], e ),
                  Concatenation( "BaumClausenInfo: failed assertion in ",
                      "inducing nonlinear representations ",
                      "(i = ", String( i ), ")\n" ) );

        fi;

        d:= Position( used, false, d );

      od;


      # Step C:
      # Compute `pilinear', `pinonlin', and `Xlist'.

      pinextlinear  := [];
      pinextnonlin1 := [];
      nextXlist1    := [];

      pinextnonlin2 := [];
      nextXlist2    := [];

      for j in [ 1 .. i-1 ] do

        pinextlinear[j]  := [];
        pinextnonlin1[j] := [];
        nextXlist1[j]    := [];

        # `cexp' describes $g_i^{g_j}$.
        cexp:= ExponentsOfPcElement( pcgs, pcgs[i]^pcgs[j], [ i .. lg ] );

        # Compute `pilinear', and the parts of `pinonlin', `Xlist'
        # arising from *linear* representations for the next step,
        # that is, compute the action of $g_j$ on `nextlinear' and
        # `nextnonlin1'.

        for k in [ 1 .. Length( linear ) ] do

          if pilinear[i][k] = k then

            # Let $F = `linear[k]'$ extend to
            # $D = D_0, D_1, \ldots, D_{p-1}$,
            # $C$ the first extension of $\pi_j(F)$.
            # We have $D( g_i^{g_j} ) = D^{g_j}(g_i) = ( C \chi^l )(g_i)$
            # where $\chi^l(g_i)$ is the $l$-th power of the chosen
            # primitive $p$-th root of unity.

            D:= nextlinear[ Dpos1[k] ];

            # `pos' is the position of $C$ in `nextlinear'.
            pos:= Dpos1[ pilinear[j][k] ];
            l:= ( (  cexp * D                   # $D( g_i^{g_j} )$
                     - nextlinear[ pos ][1] )   # $C(g_i)$
                  * p / e ) mod p;

            for u in [ 0 .. p-1 ] do
              Add( pinextlinear[j], pos + ( ( l + u * cexp[1] ) mod p ) );
            od;

          elif not IsBound( pinextnonlin1[j][ Dpos1[k] ] ) then

            # $F$ fuses with its conjugates under $g_i$,
            # the conjugating matrix describing the action of $g_j$
            # is a permutation matrix.
            # Let $D = F^{g_j}$, then the permutation corresponds to
            # the mapping between the lists
            # $[ D, (F^{g_i})^{g_j}, \ldots, (F^{g_i^{p-1}})^{g_j} ]$
            # and $[ D, D^{g_i}, \ldots, D^{g_i^{p-1}} ]$;
            # The constituents in the first list are the images of
            # the induced representation of $F$ under $g_j$,
            # and those in the second list are the constituents of the
            # induced representation of $D$.

            # While `u' runs from $1$ to $p$,
            # `pos' runs over the positions of $(F^{g_i^u})^{g_j}$ in
            # `linear'.
            # `orb' is the list of positions of the $(F^{g_j})^{g_i^u}$,
            # cyclically permuted such that the smallest entry is the
            # first.

            pinextnonlin1[j][ Dpos1[k] ]:= Dpos1[ pilinear[j][k] ];
            pos:= pilinear[j][k];
            orb:= [ pos ];
            min:= 1;
            minval:= pos;
            for u in [ 2 .. p ] do
              pos:= pilinear[i][ pos ];
              orb[u]:= pos;
              if pos < minval then
                minval:= pos;
                min:= u;
              fi;
            od;
            if 1 < min then
              orb:= Concatenation( orb{ [ min .. p ] },
                                   orb{ [ 1 .. min-1 ] } );
            fi;

            # Compute the conjugating matrix `X'.
            # Let $C$ be the stored representation $\tau_j D$
            # equivalent to $D^{g_j}$.
            # Compute the position of $C$ in `pinextnonlin1'.

            C:= nextnonlin1[ pinextnonlin1[j][ Dpos1[k] ] ];
            D:= nextnonlin1[ Dpos1[k] ];

            # `sigma' is the bijection of constituents in the restrictions
            # of $D$ and $\tau_j D$ to $G_{i-1}$.
            # More precisely, $\pi_j(\pi_i^{u-1} F) = \Phi_{\sigma(u-1)}$.
            sigma:= [];
            pos:= k;
            for u in [ 1 .. p ] do
              sigma[u]:= Position( orb, pilinear[j][ pos ] );
              pos:= pilinear[i][ pos ];
            od;

            # Compute $\pi = \sigma^{-1} (1,2,\ldots,p) \sigma$.
            pi:= [];
            pi[ sigma[p] ]:= sigma[1];
            for u in [ 1 .. p-1 ] do
              pi[ sigma[u] ]:= sigma[ u+1 ];
            od;

            # Compute the values $c_{\pi^u(0)}$, for $0 \leq u \leq p-1$.
            # Note that $c_0 = 1$.
            # (Here we encode of course the exponents.)
            constants:= [ 0 ];
            l:= 1;

            for u in [ 1 .. p-1 ] do

              # Compute $c_{\pi^u(0)}$.
              # (We have $`l' = 1 + \pi^{u-1}(0)$.)
              # Note that $B_u = [ [ 1 ] ]$ for $0\leq u\leq p-2$,
              # and $B_{p-1} = \Phi_0(g_i^p)$.

              # Next we compute the image under $A_{\pi^{u-1}(0)}$;
              # this matrix is in the $(\pi^{u-1}(0)+1)$-th column block
              # and in the $(\pi^u(0)+1)$-th row block of $D^{g_j}$.
              # Since we do not have this matrix explicitly,
              # we use the conjugate representation and the action
              # encoded by `cexp'.
              # Note the necessary initial shift because we use the
              # whole representation $D$ and not a single constituent;
              # so we shift by $\pi^u(0)+1$.
#T `perm' is nontrivial only for v = 1, this should make life easier.
              value:= 0;
              image:= pi[l];
              for v in [ 1 .. lg-i+1 ] do
                for w in [ 1 .. cexp[v] ] do
                  image:= D[v].perm[ image ];
                  value:= value + D[v].diag[ image ];
                od;
              od;

              # Next we divide by the corresponding value in
              # the image of the first standard basis vector under
              # $B_{\sigma\pi^{u-1}(0)}$.
              value:= value - C[1].diag[ sigma[l] ];
              constants[ pi[l] ]:= ( constants[l] - value ) mod e;
              l:= pi[l];

            od;

            # Put the conjugating matrix together.
            X:= rec( perm := [],
                     diag := constants );
            for u in [ 1 .. p ] do
              X.perm[ sigma[u] ]:= u;
            od;

            Assert( 2, BaumClausenInfoDebug.checkconj( pcgs, i, lg, j,
                         nextnonlin1[ Dpos1[k] ],
                         nextnonlin1[ pinextnonlin1[j][ Dpos1[k] ] ],
                         X, e ),
                  Concatenation( "BaumClausenInfo: failed assertion on ",
                      "conjugating matrices for linear repres. ",
                      "(i = ", String( i ), ")\n" ) );
            nextXlist1[j][ Dpos1[k] ]:= X;

          fi;

        od;


        # Compute the remaining parts of `pinonlin' and `Xlist' for
        # the next step, namely for those *nonlinear* representations
        # arising from *nonlinear* ones.

        nextXlist2[j]    := [];
        pinextnonlin2[j] := [];

        # `cexp' describes $g_i^{g_j}$.
        cexp:= ExponentsOfPcElement( pcgs, pcgs[i]^pcgs[j], [ i .. lg ] );

        # Compute the action of $g_j$ on `nextnonlin2'.

        for k in [ 1 .. Length( nonlin ) ] do

          if pinonlin[i][k] = k then

            # Let $F = `nonlin[k]'$ extend to
            # $D = D_0, D_1, \ldots, D_{p-1}$,
            # $C$ the first extension of $\pi_j(F)$.
            # We have $X_{j,F} \cdot F^{g_j} = \pi_j(F) \cdot X_{j,F}$,
            # thus $X_{j,F} \cdot D( g_i^{g_j} )
            # = X_{j,F} \cdot D^{g_j}(g_i)
            # = ( C \chi^l )(g_i) \cdot X_{j,F}$
            # where $\chi^l(g_i)$ is the $l$-th power of the chosen
            # primitive $p$-th root of unity.

            D:= nextnonlin2[ Dpos2[k] ];

            # `pos' is the position of $C$ in `nextnonlin2'.
            pos:= Dpos2[ pinonlin[j][k] ];

            # Find a nonzero entry in $X_{j,F} \cdot D( g_i^{g_j} )$.
            image:= Xlist[j][k].perm[1];
            value:= Xlist[j][k].diag[ image ];
            for u in [ 1 .. lg-i+1 ] do
              for v in [ 1 .. cexp[u] ] do
                image:= D[u].perm[ image ];
                value:= value + D[u].diag[ image ];
              od;
            od;

            # Subtract the corresponding value in $C(g_i) \cdot X_{j,F}$.
            C:= nextnonlin2[ pos ];
            Assert( 2, image = Xlist[j][k].perm[ C[1].perm[1] ],
                    "BaumClausenInfo: failed assertion on conj. matrices" );
            value:= value -
                ( C[1].diag[ C[1].perm[1] ] + Xlist[j][k].diag[ image ] );
            l:= ( value * p / e ) mod p;

            for u in [ 0 .. p-1 ] do
              pinextnonlin2[j][ Dpos2[k] + u ]:=
                     pos + ( ( l + u * cexp[1] ) mod p );
              nextXlist2[j][ Dpos2[k] + u ]:= Xlist[j][k];
            od;

            Assert( 2, BaumClausenInfoDebug.checkconj( pcgs, i, lg, j,
                         nextnonlin2[ Dpos2[k] ],
                         nextnonlin2[ pinextnonlin2[j][ Dpos2[k] ] ],
                         Xlist[j][k], e ),
                  Concatenation( "BaumClausenInfo: failed assertion on ",
                      "conjugating matrices for nonlinear repres. ",
                      "(i = ", String( i ), ")\n" ) );

          elif not IsBound( pinextnonlin2[j][ Dpos2[k] ] ) then

            # $F$ fuses with its conjugates under $g_i$, yielding $D$.

            dim:= Length( nonlin[k][1].diag );

            # Let $C$ be the stored representation $\tau_j D$
            # equivalent to $D^{g_j}$.
            # Compute the position of $C$ in `pinextnonlin2'.
            pinextnonlin2[j][ Dpos2[k] ]:= Dpos2[ pinonlin[j][k] ];

            C:= nextnonlin2[ pinextnonlin2[j][ Dpos2[k] ] ];
            D:= nextnonlin2[ Dpos2[k] ];

            # Compute the positions of the constituents;
            # `orb[k]' is the position of $\Phi_{k-1}$ in `nonlin'.
            pos:= pinonlin[j][k];
            orb:= [ pos ];
            min:= 1;
            minval:= pos;
            for u in [ 2 .. p ] do
              pos:= pinonlin[i][ pos ];
              orb[u]:= pos;
              if pos < minval then
                minval:= pos;
                min:= u;
              fi;
            od;
            if 1 < min then
              orb:= Concatenation( orb{ [ min .. p ] },
                                   orb{ [ 1 .. min-1 ] } );
            fi;

            # `sigma' is the bijection of constituents in the restrictions
            # of $D$ and $\tau_j D$ to $G_{i-1}$.
            # More precisely, $\pi_j(\pi_i^{u-1} F) = \Phi_{\sigma(u-1)}$.
            sigma:= [];
            pos:= k;
            for u in [ 1 .. p ] do
              sigma[u]:= Position( orb, pinonlin[j][ pos ] );
              pos:= pinonlin[i][ pos ];
            od;

            # Compute $\pi = \sigma^{-1} (1,2,\ldots,p) \sigma$.
            pi:= [];
            pi[ sigma[p] ]:= sigma[1];
            for u in [ 1 .. p-1 ] do
              pi[ sigma[u] ]:= sigma[ u+1 ];
            od;

            # Compute the positions of the constituents
            # $F_0, F_{\pi(0)}, \ldots, F_{\pi^{p-1}(0)}$.
            Forb:= [ k ];
            pos:= k;
            for u in [ 2 .. p ] do
              pos:= pinonlin[i][ pos ];
              Forb[u]:= pos;
            od;

            # Compute the values $c_{\pi^u(0)}$, for $0 \leq u \leq p-1$.
            # Note that $c_0 = 1$.
            # (Here we encode of course the exponents.)
            constants:= [ 0 ];
            l:= 1;

            for u in [ 1 .. p-1 ] do

              # Compute $c_{\pi^u(0)}$.
              # (We have $`l' = 1 + \pi^{u-1}(0)$.)
              # Note that $B_u = X_{j,\pi_j^u \Phi_0}$ for $0\leq u\leq p-2$,
              # and $B_{p-1} =
              #      \Phi_0(g_i^p) \cdot ( X_{j,\Phi_0}^{(p-1)} )^{-1}$

              # First we get the image and diagonal value of
              # the first standard basis vector under $X_{j,\pi^u(0)}$.
              image:= Xlist[j][ Forb[ pi[l] ] ].perm[1];
              value:= Xlist[j][ Forb[ pi[l] ] ].diag[ image ];

              # Next we compute the image under $A_{\pi^{u-1}(0)}$;
              # this matrix is in the $(\pi^{u-1}(0)+1)$-th column block
              # and in the $(\pi^u(0)+1)$-th row block of $D^{g_j}$.
              # Since we do not have this matrix explicitly,
              # we use the conjugate representation and the action
              # encoded by `cexp'.
              # Note the necessary initial shift because we use the
              # whole representation $D$ and not a single constituent;
              # so we shift by `dim' times $\pi^u(0)+1$.
              image:= dim * ( pi[l] - 1 ) + image;
              for v in [ 1 .. lg-i+1 ] do
                for w in [ 1 .. cexp[v] ] do
                  image:= D[v].perm[ image ];
                  value:= value + D[v].diag[ image ];
                od;
              od;

              # Next we divide by the corresponding value in
              # the image of the first standard basis vector under
              # $B_{\sigma\pi^{u-1}(0)} X_{j,\pi^{u-1}(0)}$.
              # Note that $B_v$ is in the $(v+2)$-th row block for
              # $0 \leq v \leq p-2$, in the first row block for $v = p-1$,
              # and in the $(v+1)$-th column block of $C$.
              v:= sigma[l];
              if v = p then
                image:= C[1].perm[1];
              else
                image:= C[1].perm[ v*dim + 1 ];
              fi;
              value:= value - C[1].diag[ image ];
              image:= Xlist[j][ Forb[l] ].perm[ image - ( v - 1 ) * dim ];
              value:= value - Xlist[j][ Forb[l] ].diag[ image ];
              constants[ pi[l] ]:= ( constants[l] - value ) mod e;
              l:= pi[l];

            od;

            # Put the conjugating matrix together.
            X:= rec( perm:= [],
                     diag:= [] );
            pos:= k;
            for u in [ 1 .. p ] do
              Append( X.diag, List( Xlist[j][ pos ].diag,
                                    x -> ( x + constants[u] ) mod e ) );
              X.perm{ [ ( sigma[u] - 1 )*dim+1 .. sigma[u]*dim ] }:=
                  Xlist[j][ pos ].perm + (u-1) * dim;
              pos:= pinonlin[i][ pos ];
            od;

            Assert( 2, BaumClausenInfoDebug.checkconj( pcgs, i, lg, j,
                         nextnonlin2[ Dpos2[k] ],
                         nextnonlin2[ pinextnonlin2[j][ Dpos2[k] ] ],
                         X, e ),
                  Concatenation( "BaumClausenInfo: failed assertion on ",
                      "conjugating matrices for nonlinear repres. ",
                      "(i = ", String( i ), ")\n" ) );
            nextXlist2[j][ Dpos2[k] ]:= X;

          fi;

        od;

      od;

      # Finish the update for the next index.
      BCR.linear   := nextlinear;
      BCR.pilinear := pinextlinear;

      BCR.nonlin   := Concatenation( nextnonlin1, nextnonlin2 );
      BCR.pinonlin := List( [ 1 .. i-1 ],
                       j -> Concatenation( pinextnonlin1[j],
                         pinextnonlin2[j] + Length( pinextnonlin1[j] ) ) );
      BCR.Xlist    := List( [ 1 .. i-1 ],
                    j -> Concatenation( nextXlist1[j], nextXlist2[j] ) );

      BCR.position := i-1;

      if BCR.position = 0 then
          BCR.complete := true;
      fi;

      return BCR;
end;


################################################################################
##
#F  BCLGetIrr( <BCR> )  . . . . get results of a (partical) Baum-Clausen run
##
InstallGlobalFunction( BCLGetIrr,

function(BCR)
    local
          kernel,
          linear,
          dim,
          M,
          G,
          mulmoma,        # local function  to multiply monomial matrices
          ccl,            # conjugacy classes of `G'
          tbl,            # character table of `G'
          info,           # result of `BaumClausenInfo'
          pcgs,           # value of `info.pcgs'
          lg,             # composition length
          exps,           # exponents of class representatives
          cr,             # exps sorted
          dobylevel,      # generate matrices and traces of class reps
          i, j, jj, k,    # loop variables
          t,              # intermediate representation value
          irreducibles,   # list of irreducible characters
          rep,            # loop over the representations
          l,              # list of monomial matrices
          gcd,            # g.c.d. of the exponents in `rep'
          q,              #
          o,              # order of root of unity
          e,              # exponent
          Ee,             # complex root of unity needed for `rep'
          chi,            # one character values list
          deg,            # character degree
          idmat,          # identity matrix
          perm,           # entries of monomial matrix
          diag,
          trace;          # trace of a matrix



    info := ShallowCopy(BCR);

    # Step 6: If necessary transfer the representations back to the
    #         original group.

    if  IsBound( info.hom )
       and not IsTrivial( KernelOfMultiplicativeGeneralMapping( info.hom ) ) then
      Info( InfoGroup, 2,
            "BaumClausenInfo: taking preimages in the original group" );

      kernel:= KernelOfMultiplicativeGeneralMapping( info.hom );
      k:= Pcgs( kernel );
      pcgs:= PcgsByPcSequence( ElementsFamily( FamilyObj( kernel ) ),
               Concatenation( List( pcgs,
                                    x -> PreImagesRepresentative( info.hom, x ) ),
                              k ) );
      k:= ListWithIdenticalEntries( Length( k ), 0 );

      info.linear:= List( info.linear, rep -> Concatenation( rep, k ) );

      for rep in info.nonlin do
        dim:= Length( rep[1].perm );
        M:= rec( perm:= [ 1 .. dim ],
                 diag:= [ 1 .. dim ] * 0 );
        for i in k do
          Add( rep, M );
        od;
      od;

      info.pcgs := pcgs;
      info.kernel := kernel;
    else
      info.kernel:= TrivialSubgroup( info.G );
    fi;


    mulmoma:= function( a, b )
      local prod;
      prod:= rec( perm := b.perm{ a.perm },
                  diag := [] );
      prod.diag{b.perm} := b.diag{b.perm} + a.diag;
      return prod;
    end;

    pcgs := info.pcgs{[info.position+1 .. Length(info.pcgs)]};
    G := Group(pcgs);
    pcgs:= PcgsByPcSequence( ElementsFamily( FamilyObj( G ) ), pcgs );


    tbl:= CharacterTable( G );
    ccl:= ConjugacyClasses( tbl );
    SetExponent( G, Exponent( tbl ) );
    info:=info;

    # The trivial group does not admit matrix arithmetic for evaluations.
    if IsTrivial( G ) then
      return [ Character( G, [ 1 ] ) ];
    fi;

    lg:= Length( pcgs );

    exps:= List( ccl,
                 c -> ExponentsOfPcElement( pcgs, Representative( c ) ) );

    # Compute the linear irreducibles.
    # Compute the roots of unity only once for all linear characters.
    # ($q$-th roots suffice, where $q$ divides the number of linear
    # characters and the known exponent; we do *not* compute the smallest
    # possible roots for each representation.)
    q:= Gcd( info.exponent, Length( info.linear ) );
    gcd:= info.exponent / q;
    Ee:= E(q);
    Ee:= List( [ 0 .. q-1 ], i -> Ee^i );

    irreducibles:= List( info.linear, rep ->
        Character( tbl, Ee{ ( ( exps * rep ) / gcd mod q ) + 1 } ) );

    # Compute the nonlinear irreducibles.
    if not IsEmpty( info.nonlin ) then
      cr := SortedList(exps);
      for rep in info.nonlin do
        gcd:= GcdInt( Gcd( List( rep, x -> Gcd( x.diag ) ) ), info.exponent );
        o := info.exponent / gcd;
        deg:= Length( rep[1].perm );
        idmat:= rec( perm := [ 1 .. deg ], diag := [ 1 .. deg ] * 0 );
        Add(cr[1], deg);
        l := List([1..lg], i-> idmat);
        # We go through sorted list of exponents and reuse
        # partial product from previous representative.
        for i in [2..Length(cr)] do
          j := 1;
          while cr[i-1][j] = cr[i][j] do
            j := j+1;
          od;
          for k in [cr[i-1][j]+1..cr[i][j]] do
            l[j] := mulmoma(l[j], rep[j]);
          od;
          for jj in [j+1..lg] do
            l[jj] := l[jj-1];
            for k in [1..cr[i][jj]] do
              l[jj] := mulmoma(l[jj], rep[jj]);
            od;
          od;
          # Compute the character value.
          trace:= 0*[1..o];
          perm := l[lg].perm;
          diag := l[lg].diag;
          for k in [ 1 .. deg ] do
            if perm[k] = k then
              e := (diag[k] / gcd) mod o;
              trace[e+1]:= trace[e+1] + 1;
            fi;
          od;
          # We append the character values to the exponents lists
          # and remove them (in the right order) after this loop.
          Add(cr[i], CycList(trace));
        od;
        chi := List(exps, Remove);
        Add( irreducibles, Character( tbl, chi ) );
      od;
    fi;

    # Return the result.
    return rec (
            G := G,
            irr := irreducibles,
            clreps := List(ccl, Representative)
            );
end);
