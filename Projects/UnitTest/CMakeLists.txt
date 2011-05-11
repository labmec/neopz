
set(main_SRCS Testmain.cpp)

add_executable(PZTests ${pztest} ${main_SRCS})

if(USING_BOOST)
	target_link_libraries(PZTests pthread pz boost_unit_test_framework)
else()
	target_link_libraries(PZTests pthread pz)
endif()

