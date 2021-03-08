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
        # TODOWIN32 (Gustavo 06/03/2021): this part needs testing with Visual Studio.
        # More specifically we need to check if the relative paths of the generated file tabs are
        # correct (specified here by UnitTests/${testName}).
        source_group(UnitTests/${testName} FILES ${ARGN})
    endif()
    
    target_link_libraries(${testName} PRIVATE pz Boost::unit_test_framework)
    if(USING_LOG4CXX)
      target_compile_definitions(${testName} PRIVATE LOG4CXX)
      target_link_libraries(${testName} PRIVATE ${Log4cxx_LIBRARY})
    endif()
endfunction()
