# inner product of class functions x and y
InstallGlobalFunction(IRIP,

function(GR,x,y)
	# return Sum([1..GR!.k], i ->
    #                 x[i] * ComplexConjugate(y[i]) * GR!.ccsizes[i]) / GR!.n;

     # this appears to be slightly quicker.
	 return Sum([1..GR!.k], i ->
                     x[i] * y[GR!.inverseClasses[i]] * GR!.ccsizes[i]) / GR!.n;
end
);


# inner product of class functions x and y with few entries
# one of them has all non-zero entries in classes
InstallGlobalFunction(IRIPSparse,

function(GR, classes, x,y)
    # return Sum(classes, i->
    #                 x[i] * ComplexConjugate(y[i]) * GR!.ccsizes[i]) / GR!.n;

     # this appears to be slightly quicker.
    return Sum(classes, i->
                    x[i] * y[GR!.inverseClasses[i]] * GR!.ccsizes[i]) / GR!.n;
end
);

# reduce new induced characters by the irreducibles found so far
# classes is as above
InstallGlobalFunction(IRRedSparse,
function(GR, classes)
    local mat,pos;

    if GR!.m+1 > Size(GR!.B) then
        return;
    fi;

    pos := [GR!.m+1..Size(GR!.B)];
    mat := List(pos,
                x->List(GR!.Ir, y -> IRIPSparse(GR, classes, GR!.B[x], y) ) );

    GR!.B{pos} := GR!.B{pos} - mat*GR!.Ir;
end
);

# extend the gram matrix with the scalar products of the new characters
InstallGlobalFunction(IRGramMatrixGR,

function(GR)

    local i,j,mat,b;

    b := Size(GR!.B);

    mat := NullMat(b,b);

    # keep the old inner products
    mat{[1..GR!.m]}{[1..GR!.m]} := GR!.Gram;

    for i in [GR!.m+1..b] do
        for j in [1..i] do
            mat[i][j] := IRIP(GR, GR!.B[i], GR!.B[j]);
            mat[j][i] := mat[i][j];
        od;
    od;

    GR!.Gram := mat;
end);


#############################################################################
#
# IRReduce( <GR>, ><RedTR> )
#
# reduce the induced characters by the irreducibles found so far and do LLL
# lattice reduction
#

InstallGlobalFunction(IRReduce,

function(GR, classes)

	local mat,temp,ind,I,i, Opt;

    Opt := IRGetOptions();

    IRRedSparse(GR, classes); #reduce new characters by all irreducibles
    IRGramMatrixGR(GR); # update the gram matrix

    if Opt.IRLLLOffset>0 then
        GR!.m := Size(GR!.Gram);
    else
        temp := LLLReducedGramMat(GR!.Gram, Opt.IRDelta); # LLL reduction on gram matrix
        GR!.Gram := temp.remainder;
        GR!.B := temp.transformation*GR!.B;
        temp := [];
        ind := [];
        GR!.m := Size(GR!.Gram);
        for i in [1..GR!.m] do
            if GR!.Gram[i][i]=1 then
                Add(temp,i); # find positions of characters of norm 1
            else
                Add(ind,i);
            fi;
        od;
        if Size(temp)=0 then
            Info(InfoCTUnger, 2, "Reduce: |Irr| = ", Length(GR!.Ir),
                ", dim = ", Length(GR!.Ir)+Length(GR!.Gram),
                ", det(G) = ", DeterminantMat(GR!.Gram));
            return;
        fi;
        I:=GR!.B{temp}; # irreducible characters in B (up to sign)
        if Size(ind)=0 then
            Append(GR!.Ir,I);
            GR!.Gram := [];
            Info(InfoCTUnger, 2, "Reduce: |Irr| = ", Length(GR!.Ir));
        return ;
        fi;
        for i in Reversed(temp) do
            Remove(GR!.B,i); # remove irreducible characters from B
        od;
        mat := GR!.Gram{ind}{temp};
        GR!.Gram := GR!.Gram{ind}{ind} - mat*TransposedMat(mat);
        GR!.m := Size(ind);
        GR!.B := GR!.B-mat*I;
        Append(GR!.Ir,Set(I)); # add the new irreducible characters to I
        Info(InfoCTUnger, 2, "Reduce: |Irr| = ", Length(GR!.Ir),
            ", dim = ", Length(GR!.Ir)+Length(GR!.Gram),
            ", det(G) = ", DeterminantMat(GR!.Gram));
    fi;

    return;
end);
