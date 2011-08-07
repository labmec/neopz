# @file neopz/External/Sloan/CMakeLists.txt  -- CMake file for the Sloan external module

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION include/pz)

file(GLOB sources *.cpp)

set(sloan_SRCS ${headers} ${sources} PARENT_SCOPE)
