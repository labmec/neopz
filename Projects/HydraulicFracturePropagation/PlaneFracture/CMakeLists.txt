# @file neopz/PlaneFracture/CMakeLists.txt  -- CMake file for the PlaneFracture module

file(GLOB headers *.h)
file(GLOB sources *.cpp)

if (BUILD_COMPLEX_PROJECTS)
file (GLOB to_remove_from_headers TPZPlaneFracture.h TPZJIntegral.h)
list (REMOVE_ITEM headers ${to_remove_from_headers})
file (GLOB to_remove_from_sources TPZPlaneFracture.cpp TPZJIntegral.cpp)
list (REMOVE_ITEM sources ${to_remove_from_sources})
endif()

#install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE})
set(pzplanefracture_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)
