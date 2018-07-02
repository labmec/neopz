# @file neopz/Multigrid/CMakeLists.txt  -- CMake file for the multigrid module


target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
						   $<INSTALL_INTERFACE:.> 
						  )
file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE}/Multigrid)

file(GLOB sources *.cpp)

set(pzmultigrid_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)

