if(USING_FAD)

	set(main_SRCS main.cpp)

	add_executable(ShapeFAD ${main_SRCS})

	target_link_libraries(ShapeFAD pz)

	install(TARGETS ShapeFAD ${INSTALL_TARGETS_DEFAULT_ARGS})

endif()
