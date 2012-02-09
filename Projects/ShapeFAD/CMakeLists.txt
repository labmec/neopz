if(USING_FAD)

	set(main_SRCS main.cpp)

	add_executable(ShapeFAD ${main_SRCS})

	if(USING_METIS)
		target_link_libraries(ShapeFAD pz metis)
	else()
		target_link_libraries(ShapeFAD pz)
	endif()

	install(TARGETS ShapeFAD ${INSTALL_TARGETS_DEFAULT_ARGS})

endif()
