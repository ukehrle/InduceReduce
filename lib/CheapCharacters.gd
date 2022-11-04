DeclareAttribute( "pPrimepDecompositionMaps", IsCharacterTable, "mutable" );

DeclareOperation( "pPrimepDecompositionMap", [ IsCharacterTable, IsPosInt ]  );
DeclareOperation( "pPrimepDecompositionMap", [ IsGroup, IsPosInt ]  );

DeclareGlobalFunction ("OrderCharacters");
DeclareGlobalFunction ("pPrimepDecompositionPowers");
DeclareGlobalFunction ("pPrimeCharacter");
DeclareGlobalFunction ("pPrimeRestriction");


DeclareGlobalFunction ("IRCharacterToList");
DeclareGlobalFunction ("IRListToCharacter");

DeclareGlobalFunction ("IRGenerateMoreCharacters");
DeclareGlobalFunction ("IRInitStandardChars");
