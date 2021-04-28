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
    target_link_libraries(${testName} PRIVATE test_library)
    if(PZ_LOG)
      target_link_libraries(${testName} PRIVATE ${Log4cxx_LIBRARY})
    endif()
endfunction()
