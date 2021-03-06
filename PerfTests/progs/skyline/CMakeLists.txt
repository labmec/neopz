file(GLOB source *.cpp)
file(GLOB header *.h)

add_executable(skyline-perf ${source} ${header})
target_link_libraries(skyline-perf pz)
install(TARGETS skyline-perf DESTINATION perf/bin)

