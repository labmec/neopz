# ============================================================================
#  EnableMUMPS.cmake
#  Integration of SEQUENTIAL MUMPS ([s/d/c/z]mumps) with NeoPZ
#  Based on package:
#    https://github.com/giavancini/mumps
#
# MUMPS_ROOT can be set in three ways (in order of precedence):
#   1. CMake cache variable: cmake .. -DMUMPS_ROOT=/path/to/mumps/build/local
#   2. Environment variable: export MUMPS_ROOT=/path/to/mumps/build/local
#   3. Standard system paths: /usr/local, /opt, etc. (searched by find_package)
# ============================================================================

# Set MUMPS_ROOT if not already set via cache
if(NOT MUMPS_ROOT)
    # Check environment variable
    if(DEFINED ENV{MUMPS_ROOT})
        set(MUMPS_ROOT "$ENV{MUMPS_ROOT}" CACHE PATH "MUMPS installation prefix")
        message(STATUS "[EnableMUMPS] Using MUMPS_ROOT from environment variable")
    else()
        set(MUMPS_ROOT "" CACHE PATH "MUMPS installation prefix (set via -DMUMPS_ROOT=/path or MUMPS_ROOT env var)")
    endif()
endif()

# On Apple/Clang, prefer Homebrew/MacPorts libomp when OpenMP isn't explicitly configured.
if(APPLE)
    if(NOT OpenMP_ROOT AND NOT DEFINED ENV{OpenMP_ROOT})
        if(EXISTS "/opt/homebrew/opt/libomp")
            set(OpenMP_ROOT "/opt/homebrew/opt/libomp" CACHE PATH "OpenMP (libomp) prefix" FORCE)
        elseif(EXISTS "/usr/local/opt/libomp")
            set(OpenMP_ROOT "/usr/local/opt/libomp" CACHE PATH "OpenMP (libomp) prefix" FORCE)
        elseif(EXISTS "/opt/local/lib/libomp")
            # MacPorts installs libomp in /opt/local/lib/libomp and /opt/local/include/libomp
            set(OpenMP_ROOT "/opt/local" CACHE PATH "OpenMP (libomp) prefix" FORCE)
        endif()
        if(OpenMP_ROOT)
            message(STATUS "[EnableMUMPS] OpenMP_ROOT defaulted to: ${OpenMP_ROOT}")
            list(PREPEND CMAKE_PREFIX_PATH "${OpenMP_ROOT}")
        endif()
    endif()
endif()

function(enable_mumps target)

    message(STATUS "[EnableMUMPS] Looking for MUMPS")

    # Helper to set MUMPS_DIR from MUMPS_ROOT
    if(MUMPS_ROOT)
        message(STATUS "[EnableMUMPS] MUMPS_ROOT: ${MUMPS_ROOT}")
        # Set MUMPS_DIR to help find_package locate MUMPSConfig.cmake
        # Try common subdirectory patterns
        if(EXISTS "${MUMPS_ROOT}/cmake")
            set(MUMPS_DIR "${MUMPS_ROOT}/cmake" CACHE PATH "Path to MUMPSConfig.cmake" FORCE)
        elseif(EXISTS "${MUMPS_ROOT}/lib/cmake/mumps")
            set(MUMPS_DIR "${MUMPS_ROOT}/lib/cmake/mumps" CACHE PATH "Path to MUMPSConfig.cmake" FORCE)
        elseif(EXISTS "${MUMPS_ROOT}/lib/cmake")
            set(MUMPS_DIR "${MUMPS_ROOT}/lib/cmake" CACHE PATH "Path to MUMPSConfig.cmake" FORCE)
        else()
            set(MUMPS_DIR "${MUMPS_ROOT}" CACHE PATH "Path to MUMPSConfig.cmake" FORCE)
        endif()
        if(EXISTS "${MUMPS_ROOT}/include")
            set(MUMPS_INCLUDE_DIR "${MUMPS_ROOT}/include" CACHE PATH "Path to MUMPS include directory" FORCE)
        else()
            message(WARNING "[EnableMUMPS] MUMPS include dir NOT found at: ${MUMPS_ROOT}/include")
        endif()
        message(STATUS "[EnableMUMPS] MUMPS_DIR set to: ${MUMPS_DIR}")
        message(STATUS "[EnableMUMPS] MUMPS_INCLUDE_DIR set to: ${MUMPS_INCLUDE_DIR}")

    elseif(NOT MUMPS_DIR)
        message(STATUS "[EnableMUMPS] Searching standard installation paths for MUMPS...")
        message(FATAL_ERROR
            "[EnableMUMPS] MUMPS not found in standard paths.\n"
            "Please set MUMPS_ROOT to the MUMPS installation prefix:\n\n"
            "  Option 1 - CMake command line:\n"
            "    cmake .. -DMUMPS_ROOT=/path/to/mumps/build/local\n\n"
            "  Option 2 - Environment variable (persistent):\n"
            "    export MUMPS_ROOT=/path/to/mumps/build/local\n"
            "    cmake ..\n\n"
            "  Option 3 - If MUMPS is in standard system paths:\n"
            "    Set CMAKE_PREFIX_PATH or ensure it's installed in /usr/local or /opt\n\n"
            "NOTE: Do NOT set MUMPS_DIR directly - it's an internal CMake variable.\n"
            "      Always set MUMPS_ROOT instead, which will auto-configure MUMPS_DIR.")
    else()
        message(STATUS "[EnableMUMPS] Using cached MUMPS_DIR: ${MUMPS_DIR}")
    endif()

    enable_language(C Fortran)

    # Handle OpenMP for native Apple Clang with MacPorts libomp (non-standard structure)
    # MacPorts Clang doesn't need this - it has native OpenMP support
    if(APPLE AND CMAKE_C_COMPILER_ID STREQUAL "AppleClang" AND EXISTS "/opt/local/lib/libomp")
        message(STATUS "[EnableMUMPS] Configuring OpenMP for native Apple Clang with MacPorts libomp")
        # MacPorts installs at /opt/local/lib/libomp and /opt/local/include/libomp
        if(NOT DEFINED OpenMP_C_FLAGS)
            set(OpenMP_C_FLAGS "-Xclang -fopenmp -I/opt/local/include/libomp" CACHE STRING "OpenMP C flags")
            message(STATUS "[EnableMUMPS] Setting OpenMP_C_FLAGS=\"${OpenMP_C_FLAGS}\" for Apple Clang with MacPorts libomp")
        endif()
        if(NOT DEFINED OpenMP_CXX_FLAGS)
            set(OpenMP_CXX_FLAGS "-Xclang -fopenmp -I/opt/local/include/libomp" CACHE STRING "OpenMP CXX flags")
            message(STATUS "[EnableMUMPS] Setting OpenMP_CXX_FLAGS=\"${OpenMP_CXX_FLAGS}\" for Apple Clang with MacPorts libomp")
        endif()
        if(NOT DEFINED OpenMP_C_LIB_NAMES)
            set(OpenMP_C_LIB_NAMES "omp" CACHE STRING "OpenMP library name")
            message(STATUS "[EnableMUMPS] Setting OpenMP_C_LIB_NAMES=\"${OpenMP_C_LIB_NAMES}\" for Apple Clang with MacPorts libomp")
        endif()
        if(NOT DEFINED OpenMP_CXX_LIB_NAMES)
            set(OpenMP_CXX_LIB_NAMES "omp" CACHE STRING "OpenMP library name")
            message(STATUS "[EnableMUMPS] Setting OpenMP_CXX_LIB_NAMES=\"${OpenMP_CXX_LIB_NAMES}\" for Apple Clang with MacPorts libomp")
        endif()
        if(NOT DEFINED OpenMP_omp_LIBRARY)
            set(OpenMP_omp_LIBRARY "/opt/local/lib/libomp/libomp.dylib" CACHE FILEPATH "OpenMP library path")
            message(STATUS "[EnableMUMPS] OpenMP library: ${OpenMP_omp_LIBRARY}")
        endif()
    endif()

    find_package(MUMPS CONFIG REQUIRED)
    if(NOT MUMPS_FOUND)
        message(FATAL_ERROR
            "[EnableMUMPS] MUMPS package configuration could not be loaded.\n"
            "Please verify MUMPS_ROOT or MUMPS_DIR is set correctly.")
    endif()
    if(NOT TARGET MUMPS::MUMPS)
        message(FATAL_ERROR "[EnableMUMPS] MUMPS target not found. Please check MUMPS installation.")
    endif()
    message(STATUS "[EnableMUMPS] MUMPS found: ${MUMPS_DIR}")
    message(STATUS "[EnableMUMPS] MUMPS version: ${MUMPS_VERSION}")

    # NeoPZ requires at least one of the float, double, complex and complex double ([s/d/c/z]mumps).
    message(STATUS "[EnableMUMPS] Available variants:"
        "\n  - MUMPS_s_FOUND: ${MUMPS_s_FOUND}"
        "\n  - MUMPS_d_FOUND: ${MUMPS_d_FOUND}"
        "\n  - MUMPS_c_FOUND: ${MUMPS_c_FOUND}"
        "\n  - MUMPS_z_FOUND: ${MUMPS_z_FOUND}"
    )
    if(
        NOT MUMPS_SINGLE AND NOT MUMPS_s_FOUND
        AND NOT MUMPS_DOUBLE AND NOT MUMPS_d_FOUND
        AND NOT MUMPS_COMPLEX AND NOT MUMPS_c_FOUND
        AND NOT MUMPS_COMPLEX16 AND NOT MUMPS_z_FOUND
    )
        message(FATAL_ERROR
            "[EnableMUMPS] MUMPS must be built with at least one of the following variants: [s/d/c/z]mumps")
    endif()

    # MUMPS is used to solve matrices of STATE type, so we match the MUMPS variant to STATE_TYPE,
    # not REAL_TYPE (which may differ in NeoPZ configurations where STATE != REAL).
    get_target_property(_pz_defs ${target} INTERFACE_COMPILE_DEFINITIONS)
    if(_pz_defs MATCHES "STATEfloat")
        set(_pz_state_float ON)
        message(STATUS "[EnableMUMPS] NeoPZ STATE_TYPE detected: float")
    else()
        set(_pz_state_float OFF)
        message(STATUS "[EnableMUMPS] NeoPZ STATE_TYPE detected: double (default)")
    endif()

    # Defaults: enable the real variant that matches NeoPZ STATE_TYPE
    if(_pz_state_float AND MUMPS_SINGLE)
        set(_default_single ON)
        set(_default_double OFF)
    elseif(NOT _pz_state_float AND MUMPS_DOUBLE)
        set(_default_single OFF)
        set(_default_double ON)
    else()
        set(_default_single OFF)
        set(_default_double OFF)
    endif()

    # --- User options to override defaults ---
    # Allow users to explicitly enable/disable float/double/complex/complex16 variants via cache variables.
    option(MUMPS_USE_SINGLE "Enable MUMPS float variant (smumps)" ${_default_single})
    option(MUMPS_USE_DOUBLE "Enable MUMPS double variant (dmumps)" ${_default_double})
    option(MUMPS_USE_COMPLEX "Enable MUMPS complex variant (cmumps)" OFF)
    option(MUMPS_USE_COMPLEX16 "Enable MUMPS complex16 variant (zmumps)" OFF)

    if(MUMPS_USE_SINGLE)
        message(STATUS "[EnableMUMPS] MUMPS float variant (smumps) enabled by user")
    endif()
    if(MUMPS_USE_DOUBLE)
        message(STATUS "[EnableMUMPS] MUMPS double variant (dmumps) enabled by user")
    endif()
    if(MUMPS_USE_COMPLEX)
        message(STATUS "[EnableMUMPS] MUMPS complex variant (cmumps) enabled by user")
    endif()
    if(MUMPS_USE_COMPLEX16)
        message(STATUS "[EnableMUMPS] MUMPS complex16 variant (zmumps) enabled by user")
    endif()

    # Validate: chosen variant must actually have been compiled
    if(MUMPS_USE_SINGLE AND NOT MUMPS_SINGLE)
        message(FATAL_ERROR "[EnableMUMPS] MUMPS_USE_SINGLE=ON but MUMPS single precision variant not found. Please check MUMPS installation.")
    endif()
    if(MUMPS_USE_DOUBLE AND NOT MUMPS_DOUBLE)
        message(FATAL_ERROR "[EnableMUMPS] MUMPS_USE_DOUBLE=ON but MUMPS double precision variant not found. Please check MUMPS installation.")
    endif()
    if(MUMPS_USE_COMPLEX AND NOT MUMPS_COMPLEX)
        message(FATAL_ERROR "[EnableMUMPS] MUMPS_USE_COMPLEX=ON but MUMPS complex precision variant not found. Please check MUMPS installation.")
    endif()
    if(MUMPS_USE_COMPLEX16 AND NOT MUMPS_COMPLEX16)
        message(FATAL_ERROR "[EnableMUMPS] MUMPS_USE_COMPLEX16=ON but MUMPS complex16 precision variant not found. Please check MUMPS installation.")
    endif()
    if(NOT MUMPS_USE_SINGLE AND NOT MUMPS_USE_DOUBLE AND NOT MUMPS_USE_COMPLEX AND NOT MUMPS_USE_COMPLEX16)
        message(FATAL_ERROR "[EnableMUMPS] No MUMPS variant enabled. "
            "Set at least one of MUMPS_USE_SINGLE/DOUBLE/COMPLEX/COMPLEX16=ON")
    endif()

    # ------------------------------------------------------------------------
    # 4. Compilation macros
    # ------------------------------------------------------------------------
    target_compile_definitions(${target}
        PUBLIC USING_MUMPS
        INTERFACE PZ_USING_MUMPS
    )
    if(MUMPS_USE_SINGLE)
        target_compile_definitions(${target} PUBLIC MUMPS_HAVE_SINGLE)
    endif()
    if(MUMPS_USE_DOUBLE)
        target_compile_definitions(${target} PUBLIC MUMPS_HAVE_DOUBLE)
    endif()
    if(MUMPS_USE_COMPLEX)
        target_compile_definitions(${target} PUBLIC MUMPS_HAVE_COMPLEX)
    endif()
    if(MUMPS_USE_COMPLEX16)
        target_compile_definitions(${target} PUBLIC MUMPS_HAVE_COMPLEX16)
    endif()

    if(APPLE OR NOT BUILD_SHARED_LIBS)
        target_compile_definitions(${target}
            PUBLIC MUMPS_ROOT="${MUMPS_ROOT}"
        )
        if(OpenMP_ROOT)
            target_compile_definitions(${target}
                PUBLIC OpenMP_ROOT="${OpenMP_ROOT}"
            )
        endif()
    endif()

    # ------------------------------------------------------------------------
    # 5. Include header directories and link libraries
    # INTERFACE: Propagate to downstream users but not used by this target internally

    # ------------------------------------------------------------------------
    if(NOT BUILD_SHARED_LIBS) # Still need some tests with APPLE!!
        if(gfortran) # this is to avoid the need to link with gfortran libraries manually in exterior projects
            message(STATUS "[EnableMUMPS] Linking gfortran runtime libraries for static build")
            target_link_libraries(${target} INTERFACE gfortran)
        elseif(APPLE)
            find_library(GFORTRAN_LIB gfortran PATHS /opt/homebrew/lib/gcc/current /opt/homebrew/Cellar/gcc//lib/gcc/current /opt/local/lib/gcc /opt/local/lib/libgcc)
            find_library(GOMP_LIB gomp PATHS /opt/homebrew/lib/gcc/current /opt/homebrew/Cellar/gcc//lib/gcc/current /opt/local/lib/gcc /opt/local/lib/libgcc)
            if(GFORTRAN_LIB)
                target_link_libraries(${target} PRIVATE ${GFORTRAN_LIB})
            endif()
            if(GOMP_LIB)
                target_link_libraries(${target} PRIVATE ${GOMP_LIB})
            endif()
        endif()
    endif()

    target_link_libraries(${target} PRIVATE MUMPS::MUMPS MUMPS::MPISEQ OpenMP::OpenMP_C)
    target_include_directories(${target} INTERFACE ${MUMPS_ROOT}/include)

    message(STATUS "[EnableMUMPS] MUMPS include dir: ${MUMPS_ROOT}/include")

    # Propagate MUMPS configuration variables to parent scope.
    # Used to verify if MUMPS was built with METIS support.
    set(MUMPS_metis ${MUMPS_metis} PARENT_SCOPE)

endfunction()
