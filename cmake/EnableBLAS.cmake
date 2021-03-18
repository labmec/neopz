function(enable_blas target)
    if(USING_MKL)
        set(BLA_VENDOR "Intel10_64lp") #tries to find MKL version
    endif()
    # set(BLA_VENDOR "GENERIC")
    find_package(BLAS REQUIRED)
    #target_link_libraries(${target} BLAS::BLAS) this is ok in cmake 3.18 onwards
    target_link_libraries(${target} PRIVATE ${BLAS_LIBRARIES})
    mark_as_advanced(BLAS_LIBRARIES)
    foreach(BL_LIB ${BLAS_LIBRARIES})
        if(${BL_LIB} MATCHES ".*mkl.*")
            message("Found BLAS from MKL")
            target_compile_definitions(${target} PRIVATE MKLBLAS)
            #trying to find include directory for MKL BLAS.
            get_filename_component(BL_LIB_DIR ${BL_LIB} DIRECTORY)
            if(WIN32)
              set(SEARCH_ROOT_DIR "C:")
            else()
              set(SEARCH_ROOT_DIR "/")
            endif()
            while(NOT ${BL_LIB_DIR} STREQUAL ${SEARCH_ROOT_DIR})
              find_path(BLAS_INCLUDE_DIR
                NAMES mkl_cblas.h
                PATHS ${BL_LIB_DIR}
                PATH_SUFFIXES include
                NO_DEFAULT_PATH #let us try not to mix BLAS installs
                )
              if(NOT BLAS_INCLUDE_DIR)
                #one level up...
                get_filename_component(BL_LIB_DIR_DIR ${BL_LIB_DIR} DIRECTORY)
                set(BL_LIB_DIR ${BL_LIB_DIR_DIR})
              else()
                break()
              endif()
            endwhile()
            if(NOT BLAS_INCLUDE_DIR)
              message(FATAL_ERROR "Could not find BLAS include directory")
            else()
              message(STATUS "BLAS include dirs: ${BLAS_INCLUDE_DIR}")
              target_include_directories(${target} PRIVATE ${BLAS_INCLUDE_DIR})
              mark_as_advanced(BLAS_INCLUDE_DIR)
            endif()
            break()
        endif()
    endforeach()        
endfunction()
