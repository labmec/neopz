# @file neopz/PerfUtil/CMakeLists.txt  -- CMake file for the PerfUtil module

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
						   $<INSTALL_INTERFACE:.> 
						  )

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE}/PerfUtil)

file(GLOB sources *.cpp)

set(pzperfutil_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)


