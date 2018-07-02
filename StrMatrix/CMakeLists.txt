# @file neopz/StrMatrix/CMakeLists.txt  -- CMake file for the StrMatrix module

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
						   $<INSTALL_INTERFACE:.> 
						  )

file(GLOB headers *.h)

install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE}/StrMatrix)

file(GLOB sources *.cpp)

set(pzstrmatrix_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)


