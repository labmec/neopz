
file(GLOB headers *.h)
install(FILES ${headers} DESTINATION include/pz)

file(GLOB sources *.cpp)

set(pzexternal_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)

add_subdirectory(sloan)

set(sloan_SRCS ${sloan_SRCS} PARENT_SCOPE)
