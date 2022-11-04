#############################################################################
##
#A  recipe.tst            Example package                Alexander Konovalov
##
##  To create a test file, place GAP prompts, input and output exactly as
##  they must appear in the GAP session. Do not remove lines containing 
##  START_TEST and STOP_TEST statements.
##
##  The first line starts the test. START_TEST reinitializes the caches and 
##  the global random number generator, in order to be independent of the 
##  reading order of several test files. Furthermore, the assertion level 
##  is set to 2 by START_TEST and set back to the previous value in the 
##  subsequent STOP_TEST call.
##
##  The argument of STOP_TEST may be an arbitrary identifier string.
## 
gap> START_TEST("Example package: testall.tst");

# Note that you may use comments in the test file
# and also separate parts of the test by empty lines

# Check that the global variables are defined  
gap> IsBound(DOCYCLICFIRST);
true
gap> IsBool(DOCYCLICFIRST);
true
gap> IsBound(DOCYCLICLAST);
true
gap> IsBool(DOCYCLICLAST);
true
gap> IsBound(LLLOFFSET);
true
gap> IsInt(LLLOFFSET);
true
gap> IsBound(DELTA);
true
gap> IsFloat(DELTA);
true

#############################################################################
# testing the main function
gap> CharacterTableUnger(AlternatingGroup(6));
CharacterTable( Alt( [ 1 .. 6 ] ) )

gap> CharacterTableUnger(GeneralLinearGroup(2,3));
CharacterTable( GL(2,3) )

## Each test file should finish with the call of STOP_TEST.
## The first argument of STOP_TEST should be the name of the test file.
## The second argument is redundant and is used for backwards compatibility.
gap> STOP_TEST( "testall.tst", 10000 );

#############################################################################
##
#E
