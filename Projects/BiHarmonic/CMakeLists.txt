
set(main_SRCS main.cpp)

add_executable(BiHarmonic ${main_SRCS})

if(USING_METIS)
	target_link_libraries(BiHarmonic ${KDE4_KDECORE_LIBS} pz metis pthread)
else()
	target_link_libraries(BiHarmonic ${KDE4_KDECORE_LIBS} pz pthread)
endif()

install(TARGETS BiHarmonic ${INSTALL_TARGETS_DEFAULT_ARGS})


