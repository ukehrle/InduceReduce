#
# InduceReduce: Unger's Algorithm to compute character tables of finite groups
#
# Reading the implementation part of the package.

ReadPackage( "InduceReduce", "lib/Util.gi");
ReadPackage( "InduceReduce", "lib/InduceReduce.gi");
ReadPackage( "InduceReduce", "lib/Reduce.gi");
ReadPackage( "InduceReduce", "lib/CheapCharacters.gi");

ReadPackage( "InduceReduce", "lib/profiling.gi");

if IsPackageLoaded("IO") then
    ReadPackage("InduceReduce", "lib/Parallel.gi");
fi;
