# @file neopz/Analysis/CMakeLists.txt  -- CMake file for Analysis module. 

file(GLOB headers *.h)
file(GLOB sources *.cpp)

if (BUILD_COMPLEX_PROJECTS)
file (GLOB to_remove_from_headers pzblackoilanalysis.h tpzbiharmonicestimator.h)
list (REMOVE_ITEM headers ${to_remove_from_headers})
file (GLOB to_remove_from_sources pzblackoilanalysis.cpp tpzbiharmonicestimator.cpp)
list (REMOVE_ITEM sources ${to_remove_from_sources})
endif()

install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE})
set(pzanalysis_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)
