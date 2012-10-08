########### next target ###############
file(GLOB sources *.cpp)

set(main_SRCS ${sources})

add_executable(HdivPressureTeste ${main_SRCS})

target_link_libraries(HdivPressureTeste pz)

include (CTestTestFile.cmake)