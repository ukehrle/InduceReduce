DeclareGlobalFunction( "PPart" );
DeclareGlobalFunction( "PMultiplicity" );

DeclareAttribute( "pPrimepDecompositionMaps", IsCharacterTable, "mutable" );

DeclareOperation( "pPrimepDecompositionMap", [ IsCharacterTable, IsPosInt ]  );
DeclareOperation( "pPrimepDecompositionMap", [ IsGroup, IsPosInt ]  );

DeclareGlobalFunction ("OrderCharacters");
DeclareGlobalFunction ("pPrimepDecompositionPowers");
DeclareGlobalFunction ("pPrimeCharacter");
DeclareGlobalFunction ("pPrimeRestriction");


DeclareGlobalFunction ("IRCharacterToList");
DeclareGlobalFunction ("IRListToCharacter");

DeclareGlobalFunction ("IRAddChar");
DeclareGlobalFunction ("IRAddIrrChar");
DeclareGlobalFunction ("IRGenerateMoreCharacters");
DeclareGlobalFunction ("IRInitStandardChars");
