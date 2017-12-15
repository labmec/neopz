file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE})

file(GLOB sources *.cpp)

set(pzoptimization_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)