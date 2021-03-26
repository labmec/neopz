# @file neopz/Projects/CursoPZ/GradientReconstructed/CMakeLists.txt  -- Simple project to tests gradient reconstruction by element based on least squares

file(GLOB header *.h)

file(GLOB source *.cpp)

add_executable(ApplyingGradient ${source} ${header})

target_link_libraries(ApplyingGradient pz)

