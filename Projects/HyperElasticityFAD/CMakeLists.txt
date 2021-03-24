# @file neopz/Projects/HyperElasticityFAD/CMakeLists.txt  -- CMake file for the HyperElasticityFAD project 

file(GLOB header *.h)

set(main_SRCS main.cpp)

add_executable(HyperElasticityFAD ${main_SRCS} ${header})

target_link_libraries(HyperElasticityFAD pz)
