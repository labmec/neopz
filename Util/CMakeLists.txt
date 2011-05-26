
file(GLOB headers *.h)
install(FILES ${headers} DESTINATION include/pz)

file(GLOB sources *.cpp)

set(util_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)


	set(pzsources
		${pzsources}
		${headers}
		${sources}
		PARENT_SCOPE
		)
