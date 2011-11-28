file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(main_SRCS ${sources} ${headers})

add_executable(PDElinear1D ${main_SRCS})

target_link_libraries(PDElinear1D pz pthread)

