
file(GLOB headers *.h)
file(GLOB sources *.cpp)


if(BUILD_TESTING)
	set(pztest
		${pztest}
		${headers}
		${sources}
		PARENT_SCOPE
		)
endif()

