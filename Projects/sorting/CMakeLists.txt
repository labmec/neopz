
set(sort_SRCS sort.cpp)

add_executable(sorting ${sort_SRCS})

if(USING_METIS)
	target_link_libraries(sorting ${KDE4_KDECORE_LIBS} pz metis)
else()
	target_link_libraries(sorting ${KDE4_KDECORE_LIBS} pz)
endif()

install(TARGETS sorting ${INSTALL_TARGETS_DEFAULT_ARGS})


