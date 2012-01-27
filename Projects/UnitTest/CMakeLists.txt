
set(main_SRCS Testmain.cpp)

add_executable(UnitTest_PZ ${pztest} ${main_SRCS})

if(BUILD_TESTING)
	target_link_libraries(UnitTest_PZ pthread pz boost_unit_test_framework)
else()
	target_link_libraries(UnitTest_PZ pthread pz)
endif()

