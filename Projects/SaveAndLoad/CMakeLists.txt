# @file neopz/Projects/SaveAndLoad/CMakeLists.txt  -- CMake file for save and load mesh using Write and Read methods

file(GLOB header *.h)

file(GLOB source *.cpp)

add_executable(SaveAndLoad ${source} ${header})

target_link_libraries(SaveAndLoad pz)

