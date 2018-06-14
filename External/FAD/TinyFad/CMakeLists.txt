# @file neopz/External/FAD/TinyFad/CMakeLists.txt

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
						   $<INSTALL_INTERFACE:.> 
						  )

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE}/External/FAD/TinyFad)

file(GLOB Specializations_headers Specializations/*.h)
install(FILES ${Specializations_headers} DESTINATION ${NEOPZ_INCLUDE}/External/FAD/TinyFad/Specializations)

file(GLOB sources *.cpp)

set(FAD_Tiny_SRCS ${headers} ${sources} PARENT_SCOPE)
