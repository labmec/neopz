# @file neopz/External/Sloan/CMakeLists.txt  -- CMake file for the Sloan external module

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE})

file(GLOB sources *.cpp)

set(sloan_SRCS ${headers} ${sources} PARENT_SCOPE)
