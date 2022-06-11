function(enable_omp target)
    find_package(OpenMP REQUIRED)
    include_directories(/usr/local/include)
    link_directories(/usr/local/lib)
    if (NOT OpenMP_FOUND)
        message(FATAL_ERROR "Could not find OpenMP. If OpenMP is not needed, "
                "configure the project with -DUSING_OMP=OFF")
    endif()
    #target_link_libraries(${target} PRIVATE OpenMP::OpenMP_CXX /usr/local/lib/libomp.a)
    target_link_libraries(${target} PRIVATE OpenMP::OpenMP_CXX)
    target_include_directories(${target} PRIVATE ${OpenMP_INCLUDE_DIRS})
    target_compile_definitions(${target} PRIVATE USING_OMP)
endfunction()



