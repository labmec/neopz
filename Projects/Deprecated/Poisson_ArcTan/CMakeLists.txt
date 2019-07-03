# @file neopz/Projects/Poisson_ArcTan/CMakeLists.txt  -- CMake file for using Poisson3D material bi-dimensional applying hp-adaptive

file(GLOB header *.h)

file(GLOB source *.cpp)

add_executable(Poisson_ArcTan ${source} ${header})

target_link_libraries(Poisson_ArcTan pz)

