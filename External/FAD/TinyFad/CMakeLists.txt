# @file neopz/External/FAD/TinyFad/CMakeLists.txt

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE})

file(GLOB TinyFad_headers *.h)
install(FILES ${TinyFad_headers} DESTINATION ${NEOPZ_INCLUDE}/TinyFad)

file(GLOB Specializations_headers Specializations/*.h)
install(FILES ${Specializations_headers} DESTINATION ${NEOPZ_INCLUDE}/TinyFad/Specializations)

file(GLOB sources *.cpp)

set(FAD_Tiny_SRCS ${headers} ${sources} PARENT_SCOPE)
