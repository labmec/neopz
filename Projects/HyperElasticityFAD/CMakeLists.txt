# @file neopz/Projects/HyperElasticityFAD/CMakeLists.txt  -- CMake file for the HyperElasticityFAD project 

########### next target ###############

file(GLOB header *.h)

set(main_SRCS main.cpp)

add_executable(HyperElasticity ${main_SRCS} ${header})

target_link_libraries(HyperElasticity pz pthread)

