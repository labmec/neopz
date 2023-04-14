#include <catch2/catch_session.hpp>
#include <catch2/reporters/catch_reporter_event_listener.hpp>
#include <catch2/reporters/catch_reporter_registrars.hpp>
#include <catch2/catch_test_case_info.hpp>
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

struct EventListener : Catch::EventListenerBase
{
  std::vector<std::string> failed_sections;
  std::vector<std::string> failed_but_ok_sections;
  using EventListenerBase::EventListenerBase;

  std::string lastCase="";
  bool no_fails{true};
  
  void testCaseStarting(Catch::TestCaseInfo const& testCaseInfo) final{
    lastCase = testCaseInfo.name;
    no_fails = true;
  }

  void sectionEnded(Catch::SectionStats const& sectionStats) final
  {
    if (sectionStats.assertions.failed > 0){
      /*a test case will always have one implicit section (full case).
        this logic aims to avoid printing the implicit section in case
        inner sections have already failed*/
      if(sectionStats.sectionInfo.name != lastCase || no_fails){
        failed_sections.push_back(
          lastCase
          +' '
          +sectionStats.sectionInfo.name);
        no_fails = false;
      }
    }
    else if (sectionStats.assertions.failedButOk > 0){
      /*a test case will always have one implicit section (full case).
        this logic aims to avoid printing the implicit section in case
        inner sections have already failed*/
      if(sectionStats.sectionInfo.name != lastCase || no_fails){
        failed_but_ok_sections.push_back(
          lastCase
          +' '
          +sectionStats.sectionInfo.name);
        no_fails = false;
      }
    }
  }

  void testRunEnded(Catch::TestRunStats const&) final
  {
    if(failed_but_ok_sections.size() > 0){
      std::cout<<"The following cases have failed (and were already failing before):\n";
      for(auto &s : failed_but_ok_sections){
        std::cout<<'\t'<<s<<'\n';
      }
      std::cout<<std::flush;
    }
    if(failed_sections.size() > 0){
      std::cout<<"The following cases have failed (and were NOT failing before):\n";
      for(auto &s : failed_sections){
        std::cout<<'\t'<<s<<'\n';
      }
      std::cout<<std::flush;
    }
  }
};

CATCH_REGISTER_LISTENER(EventListener)

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