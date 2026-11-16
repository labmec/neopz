# A prebuilt sequential shared MUMPS library already links its Fortran runtime.
if(TARGET MUMPS::COMMON AND NOT MUMPS_parallel AND NOT CMAKE_Fortran_COMPILER_LOADED)
    get_target_property(_mumps_library_type MUMPS::COMMON TYPE)
    if(_mumps_library_type STREQUAL "SHARED_LIBRARY")
        list(REMOVE_ITEM OpenMP_FIND_COMPONENTS Fortran)
    endif()
    unset(_mumps_library_type)
endif()

include("${CMAKE_ROOT}/Modules/FindOpenMP.cmake")
