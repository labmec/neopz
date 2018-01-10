# @file neopz/External/FAD/TinyFadET/CMakeLists.txt

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE})
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE}/TinyFadET)


file(GLOB sources *.cpp)

set(FAD_TinyET_SRCS ${headers} ${sources} PARENT_SCOPE)
