
file(GLOB headers *.h)
file(GLOB sources *.cpp)


if(BUILD_UNITTESTING)

	set(pztesttopology ${headers} ${sources} PARENT_SCOPE )

#	set(pztest
#		${pztest}
#		${headers}
#		${sources}
		
#		)

	add_executable(TestTopology ${headers} ${sources})

	target_link_libraries(TestTopology pz)

endif()

