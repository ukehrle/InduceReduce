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
##

DeclareCategory( "IsIRRecord", IsObject );

BindGlobal( "IRRecordFamily",
  NewFamily("IRRecordFamily", IsIRRecord));

BindGlobal( "IRRecordType",
  NewType(IRRecordFamily, IsIRRecord and IsAttributeStoringRep));


DeclareGlobalFunction( "CharacterTableUnger" );
DeclareGlobalFunction( "IRInit" );
DeclareGlobalFunction( "IRAddChar" );
DeclareGlobalFunction( "IRAddIrrChar" );
DeclareGlobalFunction( "IRAddStandardChars" );
DeclareGlobalFunction( "IRExcludeCovered" );
DeclareGlobalFunction( "IRExcludeCyclic" );
DeclareGlobalFunction( "IRExcludeElementary" );
DeclareGlobalFunction( "IRClassFusion" );
DeclareGlobalFunction( "IRVandermonde" );
DeclareGlobalFunction( "IRClassesCyclic" );
DeclareGlobalFunction( "IRPrintInfo" );
DeclareGlobalFunction( "IRNewElementary" );
DeclareGlobalFunction( "IRNewCyclic" );
DeclareGlobalFunction( "IRFindGroup" );
DeclareGlobalFunction( "IRHandleElementary" );
DeclareGlobalFunction( "IRIsComplete" );
DeclareGlobalFunction( "IRDoCyclics" );
DeclareGlobalFunction( "IRReduce" );
DeclareGlobalFunction( "IRtinI" );

DeclareGlobalFunction( "IRGetOpt" );
DeclareGlobalFunction( "IRGetOptions" );

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
#V IndRed
##
## a record containing subfunctions of CharacterTableUnger and InduceReduce
##
DeclareGlobalVariable( "IndRed" );

#############################################################################
##
## The info class for computations with InduceReduce
##
DeclareInfoClass("InfoCTUnger");
