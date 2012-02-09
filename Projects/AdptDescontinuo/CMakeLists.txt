
set(main_SRCS main.cpp)

add_executable(AdptDescontinuo ${main_SRCS})

if(USING_METIS)
	target_link_libraries(AdptDescontinuo ${KDE4_KDECORE_LIBS} pz metis)
else()
	target_link_libraries(AdptDescontinuo ${KDE4_KDECORE_LIBS} pz)
endif()

install(TARGETS AdptDescontinuo ${INSTALL_TARGETS_DEFAULT_ARGS})

