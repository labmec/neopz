
set(main_SRCS multicamada.cpp)

add_executable(multicamada ${main_SRCS})

if(USING_METIS)
	target_link_libraries(multicamada ${KDE4_KDECORE_LIBS} pz metis)
else()
	target_link_libraries(multicamada ${KDE4_KDECORE_LIBS} pz)
endif()

install(TARGETS multicamada ${INSTALL_TARGETS_DEFAULT_ARGS})


