# @file neopz/Refine/CMakeLists.txt  -- CMake file for the Refine module

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
						   $<INSTALL_INTERFACE:.> 
						  )

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE}/Refine)

file(GLOB sources *.cpp)

set(pzrefine_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)

# files with description of the refinement patterns

add_subdirectory(RefPatterns)
