# @file neopz/Material/Complex/CMakeLists.txt  -- CMake file for Complex Materials


target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
						   $<INSTALL_INTERFACE:.> 
						  )
file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE}/Material/Complex)

file(GLOB sources *.cpp)

set(pzmaterial_COMPLEX_SRCS ${headers} ${sources} PARENT_SCOPE)
