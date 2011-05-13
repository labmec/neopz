
file(GLOB headers *.h)
file(GLOB sources *.cpp)


if(USING_BOOST)
	set(pztest
		${pztest}
		${headers}
		${sources}
		PARENT_SCOPE
		)
endif()

