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

    IF (WIN32)
        target_link_libraries(${testName} pz ${Boost_LIBRARIES})
        # TODO (Gustavo 06/03/2021): this part needs testing with Visual Studio.
        # More specifically we need to check if the relative paths of the generated file tabs are
        # correct (specified here by UnitTests/${testName}).
        source_group(UnitTests/${testName} FILES ${ARGN})
    else()
        target_link_libraries(${testName} pz ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
    endif(WIN32)

endfunction()
