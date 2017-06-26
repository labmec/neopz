## This file should be placed in the root directory of your project.
## Then modify the CMakeLists.txt file in the root directory of your
## project to incorporate the testing dashboard.
##
## # The following are required to submit to the CDash dashboard:
##   ENABLE_TESTING()
##   INCLUDE(CTest)

set(CTEST_PROJECT_NAME "PZ")
set(CTEST_NIGHTLY_START_TIME "03:00:00 UTC")

set(CTEST_DROP_METHOD "http")
set(CTEST_DROP_SITE "labmec.org.br")
set(CTEST_DROP_LOCATION "/pz/CDash/submit.php?project=PZ")
