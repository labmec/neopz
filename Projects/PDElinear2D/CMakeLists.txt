file(GLOB headers *.h)
file(GLOB sources *.cpp)
file(GLOB sources2 ../*.cpp)

set(sources ${sources} ${sources2})

set(main_SRCS ${sources} ${headers})

add_executable(PDElinear2D ${main_SRCS})

target_link_libraries(PDElinear2D pz pthread)

