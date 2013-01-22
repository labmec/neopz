########### next target ###############
file(GLOB sources *.cpp)
file(GLOB headers *.h)



set(main_SRCS ${sources} ${headers})
add_executable(Hdiv3D ${main_SRCS})
target_link_libraries(Hdiv3D pz)

