file(GLOB header *.h)
file(GLOB source *.cpp)

add_executable(hybridmesh-perf ${source} ${header})
target_link_libraries(hybridmesh-perf pz)

