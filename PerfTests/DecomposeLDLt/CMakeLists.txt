file(GLOB source *.cpp)

file(GLOB header *.h)

add_executable(Perf-DecomposeLDLt ${source} ${header})
install(TARGETS Perf-DecomposeLDLt DESTINATION bin)
target_link_libraries(Perf-DecomposeLDLt pz)

