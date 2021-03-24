set(main_SRCS Cubo_do_Nathan.cpp)

add_executable(Cubo_do_Nathan ${main_SRCS})

target_link_libraries(Cubo_do_Nathan pz)

file(COPY cube1.txt DESTINATION .)

include (CTestTestFile.cmake)