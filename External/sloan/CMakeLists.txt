
file(GLOB headers *.h)
install(FILES ${headers} DESTINATION include/pz)


	set(pzsources
		${pzsources}
	    ${PROJECT_SOURCE_DIR}/External/sloan/diamtr.cpp
	    ${PROJECT_SOURCE_DIR}/External/sloan/gegra.cpp
	    ${PROJECT_SOURCE_DIR}/External/sloan/iorele.cpp
	    ${PROJECT_SOURCE_DIR}/External/sloan/isorti.cpp
	    ${PROJECT_SOURCE_DIR}/External/sloan/label.cpp
	    ${PROJECT_SOURCE_DIR}/External/sloan/number.cpp
	    ${PROJECT_SOURCE_DIR}/External/sloan/profi1.cpp
	    ${PROJECT_SOURCE_DIR}/External/sloan/rootls.cpp
	    ${PROJECT_SOURCE_DIR}/External/sloan/vsrtp1.cpp
		PARENT_SCOPE
		)
