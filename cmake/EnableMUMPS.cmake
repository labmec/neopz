# ============================================================================
#  EnableMUMPS.cmake
#  Integration of SEQUENTIAL MUMPS (dmumps) with NeoPZ
#  Based on package:
#    https://github.com/giavancini/mumps
# ============================================================================

function(enable_mumps target)

    message(STATUS "[EnableMUMPS] Looking for MUMPS (sequential, double)...")

    # Check if MUMPS_PREFIX is provided by user
    set(MUMPS_PREFIX "" CACHE PATH "MUMPS installation directory (root installation path)")
    if(NOT DEFINED MUMPS_PREFIX OR MUMPS_PREFIX STREQUAL "")
        message(FATAL_ERROR 
            "[EnableMUMPS] MUMPS installation directory not specified!\n"
            "Please provide:\n"
            "  -DMUMPS_PREFIX=/path/to/mumps (root installation directory)\n"
            "\n"
            "Example: cmake -DMUMPS_PREFIX=/opt/mumps ..."
        )
    endif()

    # If MUMPS_PREFIX is provided, use it to set CMAKE_PREFIX_PATH
    if(DEFINED MUMPS_PREFIX)
        list(APPEND CMAKE_PREFIX_PATH "${MUMPS_PREFIX}")
        message(STATUS "[EnableMUMPS] MUMPS_PREFIX: ${MUMPS_PREFIX}")
    endif()

    # Enable Fortran to allow OpenMP_Fortran in MUMPS targets
    enable_language(Fortran)
    
    # ------------------------------------------------------------------------
    # 1. Locate MUMPS package (requires OpenMP with Fortran)
    # ------------------------------------------------------------------------
    # Handle OpenMP on macOS with clang (similar to NeoPZ CMakeLists.txt pattern)
    if(APPLE)
        message(STATUS "[EnableMUMPS] Setting up OpenMP for macOS")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Xclang -fopenmp")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xclang -fopenmp")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")
        find_library(XCODE_OMP_LIB NAMES libiomp5.dylib omp PATHS /opt/local/lib/libomp /opt/local/lib /usr/local/lib)
        find_path(XCODE_OMP_INCLUDE omp.h PATHS /opt/local/include/libomp /opt/local/include /usr/local/include)
        if(XCODE_OMP_LIB AND XCODE_OMP_INCLUDE)
            message(STATUS "[EnableMUMPS] Found libomp: ${XCODE_OMP_LIB}")
            # Set OpenMP CMake variables so MUMPS config can find them
            set(OpenMP_C_FLAGS "-Xclang -fopenmp" CACHE STRING "C flags for OpenMP" FORCE)
            set(OpenMP_C_LIB_NAMES "omp" CACHE STRING "C OpenMP library names" FORCE)
            set(OpenMP_omp_LIBRARY "${XCODE_OMP_LIB}" CACHE FILEPATH "OpenMP library" FORCE)
            set(OpenMP_CXX_FLAGS "-Xclang -fopenmp" CACHE STRING "CXX flags for OpenMP" FORCE)
            set(OpenMP_CXX_LIB_NAMES "omp" CACHE STRING "CXX OpenMP library names" FORCE)
            set(OpenMP_Fortran_FLAGS "-fopenmp" CACHE STRING "Fortran flags for OpenMP" FORCE)
            set(OpenMP_Fortran_LIB_NAMES "omp" CACHE STRING "Fortran OpenMP library names" FORCE)
            # Set include directory
            set(OpenMP_C_INCLUDE_DIR "${XCODE_OMP_INCLUDE}" CACHE PATH "C OpenMP include directory" FORCE)
            set(OpenMP_CXX_INCLUDE_DIR "${XCODE_OMP_INCLUDE}" CACHE PATH "CXX OpenMP include directory" FORCE)
            # Mark as found
            set(OpenMP_C_FOUND TRUE CACHE BOOL "C OpenMP found" FORCE)
            set(OpenMP_CXX_FOUND TRUE CACHE BOOL "CXX OpenMP found" FORCE)
            set(OpenMP_Fortran_FOUND TRUE CACHE BOOL "Fortran OpenMP found" FORCE)
        else()
            message(FATAL_ERROR "[EnableMUMPS] Could not find libomp on macOS. Please install: sudo port install libomp")
        endif()
    else()
        find_package(OpenMP REQUIRED COMPONENTS C CXX Fortran)
    endif()
    
    find_package(MUMPS REQUIRED)

    # ------------------------------------------------------------------------
    # 2. Ensure we are NOT using parallel MUMPS
    # ------------------------------------------------------------------------
    if(MUMPS_parallel)
        message(FATAL_ERROR
            "[EnableMUMPS] MUMPS parallel detected, but this integration "
            "expects SEQUENTIAL MUMPS only."
        )
    endif()

    # ------------------------------------------------------------------------
    # 3. Ensure double precision support (dmumps)
    # ------------------------------------------------------------------------
    if(NOT MUMPS_DOUBLE)
        message(FATAL_ERROR
            "[EnableMUMPS] MUMPS was not built with DOUBLE precision (dmumps)."
        )
    endif()

    message(STATUS "MUMPS_DIR: ${MUMPS_DIR}
        parallel: ${MUMPS_parallel}
        Scotch: ${MUMPS_scotch}
        METIS: ${MUMPS_metis}
        ParMETIS: ${MUMPS_parmetis}
        OpenMP: ${MUMPS_openmp}
    ")

    message(STATUS "MUMPS_intsize64: ${MUMPS_intsize64}
        LAPACK vendor: ${MUMPS_LAPACK_VENDOR}
        Scalapack vendor: ${MUMPS_SCALAPACK_VENDOR}
        real32: ${MUMPS_SINGLE}
        real64: ${MUMPS_DOUBLE}
        complex64: ${MUMPS_COMPLEX}
        complex128: ${MUMPS_COMPLEX16}
    ")

    message(STATUS "[EnableMUMPS] Using SEQUENTIAL dmumps")
    message(STATUS "[EnableMUMPS] MUMPS version: ${MUMPS_VERSION}")

    # ------------------------------------------------------------------------
    # 4. Compilation macros
    # ------------------------------------------------------------------------
    target_compile_definitions(${target}
        PUBLIC USING_MUMPS
        INTERFACE PZ_USING_MUMPS
    )

    # ------------------------------------------------------------------------
    # 5. Include header directories
    # (MUMPS::MUMPS target already loads includes, but we keep it explicit)
    # ------------------------------------------------------------------------
    set(MUMPS_INCLUDE_DIR "${MUMPS_PREFIX}/include")
    message(STATUS "[EnableMUMPS] MUMPS include dir: ${MUMPS_INCLUDE_DIR}")
    if(DEFINED MUMPS_INCLUDE_DIR)        
        target_include_directories(${target}
            PUBLIC ${MUMPS_INCLUDE_DIR}
        )
    endif()

    # ------------------------------------------------------------------------
    # 6. Link MUMPS with OpenMP and gfortran (no parallel MPI, uses libmpiseq)
    # ------------------------------------------------------------------------
    # On macOS, we need to explicitly find and link gfortran library
    if(APPLE)
        # Add link directories for gfortran and MacPorts libraries
        target_link_directories(${target} PUBLIC /opt/local/lib /opt/local/lib/gcc14)
        
        # Find and link gfortran library explicitly
        find_library(GFORTRAN_LIB NAMES gfortran PATHS /opt/local/lib/gcc14 /opt/local/lib)
        if(NOT GFORTRAN_LIB)
            message(FATAL_ERROR "[EnableMUMPS] Could not find gfortran library. Please ensure gfortran is installed.")
        endif()
        
        # Find and link GNU OpenMP (gomp) for Fortran OpenMP support
        find_library(GOMP_LIB NAMES gomp PATHS /opt/local/lib/gcc14 /opt/local/lib)
        if(NOT GOMP_LIB)
            message(WARNING "[EnableMUMPS] Could not find libgomp. Fortran OpenMP may not work correctly.")
        endif()
        
        message(STATUS "[EnableMUMPS] Found gfortran: ${GFORTRAN_LIB}")
        if(GOMP_LIB)
            message(STATUS "[EnableMUMPS] Found gomp: ${GOMP_LIB}")
        endif()
    endif()
    
    
    target_link_libraries(${target}
        PUBLIC MUMPS::MUMPS
        "$<$<BOOL:${MPI_CXX_FOUND}>:MPI::MPI_CXX>"
        MPI::MPI_Fortran MPI::MPI_C
    )
    
    # Link gfortran and gomp on macOS, or just add -lgfortran on other platforms
    if(APPLE)
        if(GFORTRAN_LIB)
            target_link_libraries(${target} PUBLIC ${GFORTRAN_LIB})
        endif()
        if(GOMP_LIB)
            target_link_libraries(${target} PUBLIC ${GOMP_LIB})
        endif()
    else()
        target_link_libraries(${target} PUBLIC gfortran)
    endif()
    
    # Link OpenMP appropriately based on platform
    if(APPLE)
        target_link_libraries(${target} PUBLIC ${XCODE_OMP_LIB})
        target_include_directories(${target} PUBLIC ${XCODE_OMP_INCLUDE})
    else()
        target_link_libraries(${target}
            PUBLIC
            OpenMP::OpenMP_CXX
            OpenMP::OpenMP_C
            OpenMP::OpenMP_Fortran
        )
    endif()
    
    # always link MPI (MPI::MPI_Fortran and MPI::MPI_C) as we use the libmpiseq if not MUMPS_parallel
    
    # Force OpenMP flags in compilation (use appropriate syntax for compiler)
    if(APPLE)
        target_compile_options(${target} PUBLIC -Xclang -fopenmp)
        target_link_options(${target} PUBLIC -Xclang -fopenmp)
    else()
        target_compile_options(${target} PUBLIC -fopenmp)
        target_link_options(${target} PUBLIC -fopenmp)
    endif()

    message(STATUS "[EnableMUMPS] MUMPS sequential enabled for target: ${target}")

endfunction()
