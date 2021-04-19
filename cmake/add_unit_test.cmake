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
    target_compile_definitions(${testName} PRIVATE CATCH_CONFIG_MAIN)
    target_link_libraries(${testName} PRIVATE pz Catch2::Catch2)
    if(PZ_LOG)
      target_link_libraries(${testName} PRIVATE ${Log4cxx_LIBRARY})
    endif()
    if (WIN32)
        # TODOWIN32 (Gustavo 06/03/2021): this part needs testing with Visual Studio.
        # More specifically we need to check if the relative paths of the generated file tabs are
        # correct (specified here by UnitTests/${testName}).
        source_group(UnitTests/${testName} FILES ${ARGN})
        foreach(file ${CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS})
            get_filename_component(fileName ${file} NAME)
            add_custom_command(TARGET ${testName} POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different
                ${file}
                ${CMAKE_CFG_INTDIR}/${fileName})
        endforeach()
        add_custom_command(TARGET ${testName} POST_BUILD      COMMAND ${CMAKE_COMMAND} -E copy_if_different      "${PROJECT_BINARY_DIR}/$<CONFIG>/pz.dll"      $<TARGET_FILE_DIR:${testName}>/pz.dll)
    endif()
endfunction()
