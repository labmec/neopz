file(GLOB source *.cpp)

file(GLOB header *.h)

add_executable(Perf-DecomposeLDLt ${source} ${header})

target_link_libraries(Perf-DecomposeLDLt pz pthread)

