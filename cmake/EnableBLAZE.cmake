function(enable_blaze target)
    if (NOT USING_MKL AND NOT USING_LAPACK)
        message(FATAL_ERROR "Blaze library requires either USING_MKL or USING_LAPACK enabled!")
    endif()
    find_package(blaze REQUIRED)
    target_compile_definitions(${target} PRIVATE USING_BLAZE)
    target_link_libraries(${target} PRIVATE blaze::blaze)
    # Blaze docs enforce this line. Mandatory under Windows. In future, should include "INTERFACE" modifier.
endfunction()
