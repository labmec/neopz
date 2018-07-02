# Project to apply hp-adaptive process to Poisson3D material on hexahedra with a shock problem

file(GLOB header *.h ../LibRefine/*.h)
file(GLOB source *.cpp ../LibRefine/*.cpp)

set(HPAdaptive_SRCS ${header} ${source})

add_executable(Poisson3D_Shock ${HPAdaptive_SRCS})

target_link_libraries(Poisson3D_Shock pz)

