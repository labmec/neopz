
file(GLOB headers *.h)
install(FILES ${headers} DESTINATION include/pz)

set(pzsources
	${pzsources}
	${headers}
	PARENT_SCOPE
	)

