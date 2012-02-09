

set(main_SRCS main.cpp)

add_executable(BRUReader ${main_SRCS})

if(USING_METIS)
	target_link_libraries(BRUReader pz metis)
else()
	target_link_libraries(BRUReader pz)
endif()

install(TARGETS BRUReader ${INSTALL_TARGETS_DEFAULT_ARGS})


