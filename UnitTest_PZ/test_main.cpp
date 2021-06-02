#define CATCH_CONFIG_RUNNER
#include <catch2/catch.hpp>
#ifdef PZ_LOG
#include "pzlog.h"
#endif

int main( int argc, char* argv[] ) {
  // global setup...
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
  int result = Catch::Session().run( argc, argv );

  // global clean-up...

  return result;
}

/*
TEMPLATE_TEST_CASE_SYNTAX:

"name",
"[my-tags][as-many-tags-as-i-want]",
type1,
type2,
...,
typen

Test name must be unique, tags can be repeated. 
Then you can run just one test as, e.g.:
./TestStruct -n "Assemble known matrix"


At the moment of writing, the tags are not organised.

I propose we use tags as:

TEMPLATE_TEST_CASE("Compare parallel and serial matrices","[struct_tests][struct][multithread]",
                   TPZStructMatrixOR<STATE>,TPZStructMatrixOT<STATE>)

first tag: test name
other tags: classes tested (struct, analysis, etc), important aspects (multithread).

Check existing tags before adding yours!
*/