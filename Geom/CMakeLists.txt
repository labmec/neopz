
file(GLOB headers *.h)
install(FILES ${headers} DESTINATION include/pz)

file(GLOB sources *.cpp)

set(pzgeom_STAT_SRCS ${sources} ${headers} PARENT_SCOPE)

	set(pzsources
		${pzsources}
		${headers}
		${sources}
		PARENT_SCOPE
		)

