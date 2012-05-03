# @file neopz/Pre/CMakeLists.txt  -- CMake file for the Pre module

file(GLOB headers *.h)

if (BUILD_COMPLEX_PROJECTS)
file (GLOB to_remove_from_headers pzreadmeshhr.h)
list (REMOVE_ITEM headers ${to_remove_from_headers})
endif()

install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE})

file(GLOB sources *.cpp)

if (BUILD_COMPLEX_PROJECTS)
file (GLOB to_remove_from_sources pzreadmeshhr.cpp)
list (REMOVE_ITEM sources ${to_remove_from_sources})
endif()

set(pzpre_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)

