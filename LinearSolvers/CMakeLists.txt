# @file neopz/LinearSolvers/CMakeLists.txt  -- CMake file for the LinearSolvers module

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
						   $<INSTALL_INTERFACE:.> 
						  )

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE}/LinearSolvers)

file(GLOB sources *.cpp)

set(pzlinearsolvers_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)


