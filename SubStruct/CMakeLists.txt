# @file neopz/SubStruct/CMakeLists.txt  -- CMake file for the SubStruct module

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
						   $<INSTALL_INTERFACE:.> 
						  )

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE}/SubStruct)

file(GLOB sources *.cpp)

set(pzsubstruct_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)


