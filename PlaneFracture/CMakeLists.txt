# @file neopz/PlaneFracture/CMakeLists.txt  -- CMake file for the PlaneFracture module

file(GLOB headers *.h)

if (BUILD_COMPLEX_PROJECTS)
file (GLOB to_remove_from_headers TPZPlaneFracture.h)
list (REMOVE_ITEM headers ${to_remove_from_headers})
endif()

install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE})

file(GLOB sources *.cpp)

if (BUILD_COMPLEX_PROJECTS)
file (GLOB to_remove_from_sources TPZPlaneFracture.cpp)
list (REMOVE_ITEM sources ${to_remove_from_sources})
endif()

set(pzplanefracture_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)


