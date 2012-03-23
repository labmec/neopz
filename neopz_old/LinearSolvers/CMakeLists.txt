# @file neopz/LinearSolvers/CMakeLists.txt  -- CMake file for the LinearSolvers module

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE})

file(GLOB sources *.cpp)

set(pzlinearsolvers_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)


