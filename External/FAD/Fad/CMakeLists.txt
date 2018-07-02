# @file neopz/External/FAD/Fad/CMakeLists.txt

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
						   $<INSTALL_INTERFACE:.> 
						  )

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE}/External/FAD/Fad)

file(GLOB sources *.cpp)

set(FAD_fad_SRCS ${headers} ${sources} PARENT_SCOPE)
