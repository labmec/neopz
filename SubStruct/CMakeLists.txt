
file(GLOB headers *.h)
install(FILES ${headers} DESTINATION include/pz)

file(GLOB sources *.cpp)

set(pzsources
	${pzsources}
	${headers}
	${sources}
	PARENT_SCOPE
)
set(pzsubstruct_STAT_SRCS ${sources} ${headers} PARENT_SCOPE)
