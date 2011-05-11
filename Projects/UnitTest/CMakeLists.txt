
set(main_SRCS Testmain.cpp)

add_executable(PZTests ${pztest} ${main_SRCS})

target_link_libraries(PZTests pthread pz boost_unit_test_framework)

