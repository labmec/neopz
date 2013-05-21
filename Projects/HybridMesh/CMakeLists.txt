
file(GLOB header *.h)

file(GLOB source *.cpp)

add_executable(HybridMesh ${source} ${header})

target_link_libraries(HybridMesh pz)

