file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(main_SRCS ${sources} ${headers})

add_executable(CursoPZ_Laplacian2D ${main_SRCS})

target_link_libraries(CursoPZ_Laplacian2D pz)
