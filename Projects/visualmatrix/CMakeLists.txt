# @file neopz/Projects/visualmatrix/CMakeLists.txt  -- CMake file for the VisualMatrix tutorial 

file(GLOB source *.cpp)
file(GLOB header *.h)

set(visualmatrix ${visualmatrix} ${source} ${header} PARENT_SCOPE)

add_executable(VisualMatrix ${source} ${header})

target_link_libraries(VisualMatrix pz pthread)
