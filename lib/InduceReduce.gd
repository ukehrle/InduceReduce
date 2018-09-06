#
# InduceReduce.gd        The InduceReduce package         Jonathan Gruber
#
# Declarations
#
# Copyright (C) 2018     Jonathan Gruber
#
#############################################################################
##
#F CharacterTableUnger( <G> )
## 
## Computes the character table of a finite group using Unger's algorithm
##
DeclareGlobalFunction( "CharacterTableUnger" );

#############################################################################
##
#F IndRedInit( <G> )
##
## Initialization of the elgorithm returns a record which contains all
## relevant data concerning the group G
##
DeclareGlobalFunction( "IndRedInit" );

#############################################################################
##
#F IndRedGroupTools( )
##
## returns a record with different functions used for the algorithm
##
DeclareGlobalFunction( "IndRedGroupTools" );

#############################################################################
##
#F IndRedReduceTools( )
##
## returns a record with different functions required for lattice reduction
##
DeclareGlobalFunction( "IndRedReduceTools" );

#############################################################################
##
#F IndRedFindElementary( <GR> )
##
## finds the next elementary subgroup and puts it in the component Elementary
## of the record GR.
##
DeclareGlobalFunction( "IndRedFindElementary" );

#############################################################################
##
#F IndRedInitElementary( <GR> , <TR> )
##
## initialize data concerning the elementary subgroups and exclude its
## subgroups and their conjugates from the computations yet to come.
##
DeclareGlobalFunction( "IndRedInitElementary" );

#############################################################################
#
# InduceReduce.gd        The InduceReduce package         Jonathan Gruber
#
# Declarations
#
# Copyright (C) 2018     Jonathan Gruber
#
#############################################################################
##
#F CharacterTableUnger( <G> [, <Options>] )
## 
## Computes the character table of a finite group using Unger's algorithm
##
DeclareGlobalFunction( "CharacterTableUnger" );

#############################################################################
##
#F IndRedInit( <G> )
##
## Initialization of the elgorithm returns a record which contains all
## relevant data concerning the group G
##
DeclareGlobalFunction( "IndRedInit" );

#############################################################################
##
#F IndRedGroupTools( )
##
## returns a record with different functions used for the algorithm
##
DeclareGlobalFunction( "IndRedGroupTools" );

#############################################################################
##
#F IndRedReduceTools( )
##
## returns a record with different functions required for lattice reduction
##
DeclareGlobalFunction( "IndRedReduceTools" );

#############################################################################
##
#F IndRedFindElementary( <GR> , <Opt> )
##
## finds the next elementary subgroup and puts it in the component Elementary
## of the record GR.
##
DeclareGlobalFunction( "IndRedFindElementary" );

#############################################################################
##
#F IndRedInitElementary( <GR> , <TR> )
##
## initialize data concerning the elementary subgroups and exclude its
## subgroups and their conjugates from the computations yet to come.
##
DeclareGlobalFunction( "IndRedInitElementary" );

#############################################################################
##
#F IndRedInduceCyc( <GR> )
##
## induce characters from all cyclic subgroups
##
DeclareGlobalFunction( "IndRedInduceCyc" );

#############################################################################
##
#F IndRedInduce( <GR> )
##
## Induce all irreducible characters of the elementary subgroup in 
## GR.Elementary to the group GR.G and add them to GR.B
##
DeclareGlobalFunction( "IndRedInduce" );

#############################################################################
#
# IndRedReduce( <GR>, ><RedTR> )
#
# reduce the induced characters by the irreducibles found so far and do LLL
# lattice reduction
#
DeclareGlobalFunction( "IndRedReduce" );

#############################################################################
##
#F InduceReduce( <GR> , <Opt> )
##
## computes the irreducible characters of the group GR.G (upto sign) and puts
## them in the component GR.Ir of GR, returns GR.Ir
##
DeclareGlobalFunction( "InduceReduce" );

#############################################################################
##
#F IndRedtinI( <GR> )
##
## adjusts the signs of the characters and permutes them to the right
## ordering of conjugacy classes
##
DeclareGlobalFunction( "IndRedtinI" );

#############################################################################
##
## The info class for computations with InduceReduce
##
DeclareInfoClass("InfoCTUnger");
