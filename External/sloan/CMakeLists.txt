
file(GLOB headers *.h)
install(FILES ${headers} DESTINATION include/pz)

file(GLOB sources *.cpp)

set(sloan_SRCS ${headers} ${sources} PARENT_SCOPE)
