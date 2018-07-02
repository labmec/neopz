# @file neopz/Shape/CMakeLists.txt  -- CMake file for the Shape module

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
						   $<INSTALL_INTERFACE:.> 
						  )

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE}/Shape)

file(GLOB sources *.cpp)

set(pzshape_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)




