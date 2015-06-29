file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(main_SRCS ${sources} ${headers})

add_executable(TesteFranComplexo2D ${main_SRCS})

target_link_libraries(TesteFranComplexo2D pz)