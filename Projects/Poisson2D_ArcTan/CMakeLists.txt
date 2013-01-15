# @file neopz/Projects/Poisson2D_ArcTan/CMakeLists.txt  -- CMake file for using Poisson3D material bi-dimensional applying hp-adaptive

file(GLOB header *.h)

file(GLOB source *.cpp)

add_executable(Poisson2D_ArcTan ${source} ${header})

target_link_libraries(Poisson2D_ArcTan pz)

