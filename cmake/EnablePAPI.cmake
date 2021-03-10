function(enable_papi target)
    find_package(PAPI REQUIRED)
    target_link_libraries(${target} PRIVATE ${PAPI_LIB})
endfunction()
