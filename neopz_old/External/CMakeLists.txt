# @file neopz/External/CMakeLists.txt  -- CMake file for the External module

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE})

file(GLOB sources *.cpp)

set(pzexternal_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)

add_subdirectory(sloan)

set(sloan_SRCS ${sloan_SRCS} PARENT_SCOPE)
