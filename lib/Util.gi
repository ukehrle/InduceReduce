InstallGlobalFunction( PPart,

function(n,q)
    local res;
    res:=1;
    while RemInt(n,res)=0 do
        res:=res*q;
    od;
    return res/q;
end);

InstallGlobalFunction( PMultiplicity,

function(n,p)
    local res;
    res := 0;
    while IsInt(n/p) do
        n := n/p;
        res := res +1;
    od;
    return res;
end);
