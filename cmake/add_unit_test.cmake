#
# Create a Boost unit test with a respective CTest test named 'testName'.
# All remaining (unnamed) arguments are treated as the source files for the test case,
# accessed by the variable 'ARGN'.
#
# Usage:
#     add_unit_test(MyTest source.cpp header.h ... )
#
function(add_unit_test testName)

    add_test(${testName} ${testName})
    add_executable(${testName} ${ARGN})
    target_link_libraries(${testName} PUBLIC test_library)
    # not working properly, log4cxx calls are detected as tests
    # include(${Catch2_SOURCE_DIR}/extras/Catch.cmake)
    # catch_discover_tests(${testName})
endfunction()
