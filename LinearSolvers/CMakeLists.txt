
file(GLOB headers *.h)
install(FILES ${headers} DESTINATION include/pz)

set(pzsources
	${pzsources}
	${headers}
	PARENT_SCOPE
	)

	set(pzlinearsolvers_STAT_SRCS ${sources} ${headers} PARENT_SCOPE)
