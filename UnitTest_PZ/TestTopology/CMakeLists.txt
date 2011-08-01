
file(GLOB headers *.h)
file(GLOB sources *.cpp)


if(BUILD_UNITTESTING)

	set(pztestlocal ${headers} ${sources} )

	set(pztest
		${pztest}
		${headers}
		${sources}
		PARENT_SCOPE
		)

	add_executable(TestTopology ${pztestlocal})

	target_link_libraries(TestTopology pz)

endif()

