# @file neopz/External/FAD/TinyFad/CMakeLists.txt

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE})

file(GLOB sources *.cpp)

set(FAD_Tiny_SRCS ${headers} ${sources} PARENT_SCOPE)
