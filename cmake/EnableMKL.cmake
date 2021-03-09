function(enable_mkl target)
    if (WIN32)
        set(MKL_ROOT "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows")
        find_path(MKL_INCLUDE NAMES mkl.h PATHS ${MKL_ROOT}/mkl/include)
        find_library(MKL_LIB_INTEL NAMES mkl_rt.lib PATHS ${MKL_ROOT}/mkl/lib/intel64_win)
        find_library(MKL_LIB_CORE NAMES mkl_core.lib PATHS ${MKL_ROOT}/mkl/lib/intel64_win)
        find_library(MKL_LIB_THREAD NAMES mkl_intel_thread.lib PATHS ${MKL_ROOT}/mkl/lib/intel64_win)
        find_library(COMPOSER_OMP NAMES libiomp5md.lib PATHS ${MKL_ROOT}/compiler/lib/intel64_win)
    else()
        find_path(MKL_INCLUDE NAMES mkl.h PATHS ${SEARCH_DIRS} /opt/intel/mkl/include /softwares/intel/mkl/include)
        find_library(MKL_LIB_INTEL NAMES libmkl_intel_lp64.so libmkl_intel_lp64.dylib PATHS
                     ${SEARCH_DIRS} /opt/intel/mkl/lib /opt/intel/mkl/lib/intel64/ /softwares/intel/mkl/lib /softwares/intel/mkl/lib/intel64)
        find_library(MKL_LIB_CORE NAMES libmkl_core.so libmkl_core.dylib PATHS ${SEARCH_DIRS}
                     /opt/intel/mkl/lib /opt/intel/mkl/lib/intel64/ /softwares/intel/mkl/lib /softwares/intel/mkl/lib/intel64/)
        find_library(MKL_LIB_THREAD NAMES libmkl_intel_thread.so libmkl_intel_thread.dylib PATHS
                     ${SEARCH_DIRS} /opt/intel/mkl/lib /opt/intel/mkl/lib/intel64/ /softwares/intel/mkl/lib /softwares/intel/mkl/lib/intel64/)
        find_library(COMPOSER_OMP NAMES libiomp5.so libiomp5.dylib PATHS ${SEARCH_DIRS}
                     /opt/intel/composer_xe/compiler/lib
                     /opt/intel/composerxe/lib/intel64
                     /opt/intel/lib /softwares/intel/lib/intel64 /opt/intel/mkl/lib/intel64/
                     /opt/intel/compilers_and_libraries/linux/lib/intel64_lin )
    endif()

    if(MKL_INCLUDE-NOTFOUND)
        set (MKL_INCLUDE "" CACHE PATH "Directory where mkl.h can be found")
    else()
        target_include_directories(${target} PRIVATE ${MKL_INCLUDE})
    endif()

    if(MKL_LIB-NOTFOUND)
        set (MKL_LIB_INTEL "" CACHE PATH "Directory where the mkl library can be found")
    else()
        target_link_libraries(${target} PRIVATE ${MKL_LIB_INTEL})
        target_link_libraries(${target} PRIVATE ${MKL_LIB_CORE})
        target_link_libraries(${target} PRIVATE ${MKL_LIB_THREAD})
        target_link_libraries(${target} PRIVATE ${COMPOSER_OMP})
        target_compile_definitions(${target} PRIVATE MKLBLAS)
        target_compile_definitions(${target} PRIVATE MKLLAPACK)
    endif()
endfunction()
