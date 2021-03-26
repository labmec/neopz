# @file neopz/Publications/CMakeLists.txt  -- CMake file for the Publications module

if (NOT BUILD_COMPLEX_PROJECTS) ##projetos nao compativeis com state=complex<double>
    add_pz_internal_target(NAME hdiv2d SOURCES hdiv2dpaper201504.cpp hdiv2dpaper201504.h)
    add_pz_internal_target(NAME hdiv3d SOURCES hdiv3dpaper201504.cpp hdiv3dpaper201504.h)
    add_pz_internal_target(NAME hdivCurved SOURCES hdivCurvedJCompAppMath.cpp hdivCurvedJCompAppMath.h)
endif ()
