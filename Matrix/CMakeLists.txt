
file(GLOB headers *.h)
install(FILES ${headers} DESTINATION include/pz)

file(GLOB sources *.cpp)


	set(pzsources
		${pzsources}
		${headers}
		${sources}
		PARENT_SCOPE
		)

if(BUILD_TESTING)
	add_subdirectory(Test)

	set(pztest
		${pztest}
		PARENT_SCOPE
		)
endif()