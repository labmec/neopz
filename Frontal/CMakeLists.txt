# @file neopz/Frontal/CMakeLists.txt  -- CMake file for the Frontal module

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
						   $<INSTALL_INTERFACE:.> 
						  )

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE}/Frontal)

file(GLOB sources *.cpp)

set(pzfrontal_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)

