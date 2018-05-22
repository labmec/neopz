# @file neopz/SpecialMaps/CMakeLists.txt  -- CMake file for the SpecialMaps module

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
						   $<INSTALL_INTERFACE:.> 
						  )

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE}/SpecialMaps)

file(GLOB sources *.cpp)

set(pzspecialmaps_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)



