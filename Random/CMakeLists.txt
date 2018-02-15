file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE})

file(GLOB sources *.cpp)

set(pzrandom_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)