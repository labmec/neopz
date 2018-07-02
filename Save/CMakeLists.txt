# @file neopz/Save/CMakeLists.txt  -- CMake file for the Save module

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
						   $<INSTALL_INTERFACE:.> 
						  )

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE}/Save)

file(GLOB sources *.cpp)

set(pzsave_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)


