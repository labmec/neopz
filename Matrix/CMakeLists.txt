# @file neopz/Matrix/CMakeLists.txt  -- CMake file for the Matrix module

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
						   $<INSTALL_INTERFACE:.> 
						  )

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE}/Matrix)

file(GLOB sources *.cpp)

set(pzmatrix_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)

