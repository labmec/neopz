# @file neopz/Projects/Poisson2D_ArcTan_Disc/CMakeLists.txt  -- CMake file for using Poisson3D material bi-dimensional applying hp-adaptive

file(GLOB header *.h)

file(GLOB source *.cpp)

add_executable(Poisson2D_ArcTan_Disc ${source} ${header})

target_link_libraries(Poisson2D_ArcTan_Disc pz)

