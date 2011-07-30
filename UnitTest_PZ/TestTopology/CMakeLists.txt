
file(GLOB headers *.h)
file(GLOB sources *.cpp)


if(BUILD_UNITTESTING)
	set(pztest
		${pztest}
		${headers}
		${sources}
		)

	add_executable(TestTopology ${pztest})

	target_link_libraries(TestTopology pz)

endif()

