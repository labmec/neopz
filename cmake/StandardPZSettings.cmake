# Set a default build type if none was specified
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to 'RelWithDebInfo' as none was specified.")
    set(CMAKE_BUILD_TYPE
        RelWithDebInfo
        CACHE STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui, ccmake
    set_property(
        CACHE CMAKE_BUILD_TYPE
        PROPERTY STRINGS
        "Debug"
        "Release"
        "MinSizeRel"
        "RelWithDebInfo")
endif ()

#To select type of REAL
if (NOT REAL_TYPE)
    message(STATUS "Setting REAL type to 'double' as none was specified.")
    set(REAL_TYPE
        double
        CACHE STRING "Choose the type of REAL numbers." FORCE)
    # Set the possible values of build type for cmake-gui, ccmake
    set_property(
        CACHE REAL_TYPE
        PROPERTY STRINGS
        "double"
        "float"
        "long double"
        "pzfpcounter")
endif ()

if (REAL_TYPE STREQUAL "double")
    set(REAL_TYPE_DEF "REALdouble")
elseif (REAL_TYPE STREQUAL "float")
    set(REAL_TYPE_DEF "REALfloat")
elseif (REAL_TYPE STREQUAL "long double")
    set(REAL_TYPE_DEF "REALlongdouble")
elseif (REAL_TYPE STREQUAL "pzfpcounter")
    set(REAL_TYPE_DEF "REALpzfpcounter")
else()
    message(FATAL_ERROR "Please specify a valid type for REAL from the options: "
            "'float', 'double', 'long double', 'pzfpcounter'.")
endif()

#To select type of STATE
if (NOT STATE_TYPE)
    message(STATUS "Setting STATE type to 'double' as none was specified.")
    set(STATE_TYPE
        double
        CACHE STRING "Choose the type of the STATE variable." FORCE)
    # Set the possible values of build type for cmake-gui, ccmake
    set_property(
        CACHE STATE_TYPE
        PROPERTY STRINGS
        "double"
        "float"
        "long double")
endif ()

if (STATE_TYPE STREQUAL "double")
    set(STATE_TYPE_DEF "STATEdouble")
elseif (STATE_TYPE STREQUAL "float")
    set(STATE_TYPE_DEF "STATEfloat")
elseif (STATE_TYPE STREQUAL "long double")
    set(STATE_TYPE_DEF "STATElongdouble")
else()
    message(FATAL_ERROR "Please specify a valid type for STATE from the options: "
            "'double', 'float', 'long double'")
endif()

message(STATUS "NeoPZ configuration:"
        "\n - Build type: " ${CMAKE_BUILD_TYPE}
        "\n - REAL type: " ${REAL_TYPE}
        "\n - STATE type: " ${STATE_TYPE} "\n")

#specifying the path to neopz source code
set(PZSOURCEDIR ${PROJECT_SOURCE_DIR})
#specify where the refinement patterns can be found (for build)
set(PZ_REFPATTERN_DIR ${PROJECT_SOURCE_DIR}/Refine/RefPatterns)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
