
set(main_SRCS main.cpp)

add_executable(AltoDes ${main_SRCS})

if(USING_METIS)
	target_link_libraries(AltoDes ${KDE4_KDECORE_LIBS} pz metis)
else()
	target_link_libraries(AltoDes ${KDE4_KDECORE_LIBS} pz)
endif()

install(TARGETS AltoDes ${INSTALL_TARGETS_DEFAULT_ARGS})

