#
# InduceReduce: Unger's Algorithm to compute character tables of finite groups
#
# Reading the declaration part of the package.
#

ReadPackage( "InduceReduce", "lib/Util.gd");
ReadPackage( "InduceReduce", "lib/InduceReduce.gd");
ReadPackage( "InduceReduce", "lib/Reduce.gd");
ReadPackage( "InduceReduce", "lib/CheapCharacters.gd");

if IsPackageLoaded("IO") then
    ReadPackage("InduceReduce", "lib/Parallel.gd");
fi;
