# @file neopz/External/FAD/Fad/CMakeLists.txt

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE})

file(GLOB sources *.cpp)

set(FAD_fad_SRCS ${headers} ${sources} PARENT_SCOPE)
