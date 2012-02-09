
set(RefPat_SRCS main.cpp)

add_executable(Test_RefPattern ${RefPat_SRCS})

if(USING_METIS)
	target_link_libraries(Test_RefPattern pz metis)
else()
	target_link_libraries(Test_RefPattern pz)
endif()

install(TARGETS Test_RefPattern ${INSTALL_TARGETS_DEFAULT_ARGS})

