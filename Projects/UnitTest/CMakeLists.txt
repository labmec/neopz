
set(main_SRCS Testmain.cpp)

add_executable(PZTests ${pztest} ${main_SRCS})

if(BUILD_TESTING)
	target_link_libraries(PZTests pthread pz boost_unit_test_framework)
else()
	target_link_libraries(PZTests pthread pz)
endif()

