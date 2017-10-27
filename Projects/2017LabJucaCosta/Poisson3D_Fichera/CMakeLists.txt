# Project to apply hp-adaptive process to Poisson3D material on Fichera corner

file(GLOB header *.h)
file(GLOB source *.cpp)

set(HPAdaptive_SRCS ${header} ${source})

add_executable(Poisson3D_Fichera ${HPAdaptive_SRCS})

target_link_libraries(Poisson3D_Fichera pz)

