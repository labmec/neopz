# @file neopz/External/Sloan/CMakeLists.txt  -- CMake file for the Sloan external module

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
						   $<INSTALL_INTERFACE:.> 
						  )

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE}/External/Sloan)

file(GLOB sources *.cpp)

set(sloan_SRCS ${headers} ${sources} PARENT_SCOPE)
