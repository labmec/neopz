
file(GLOB headers *.h)
install(FILES ${headers} DESTINATION include)

file(GLOB sources *.cpp)

########### next target ###############

set(swxLib_STAT_SRCS
	${headers}
	${sources}
	${swxLib_STAT_SRCS}
	PARENT_SCOPE
	)

