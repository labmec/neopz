file(GLOB source *.cpp)
file(GLOB header *.h)

add_executable(skylmat-perf ${source} ${header})
target_link_libraries(skylmat-perf pz)
install(TARGETS skylmat-perf DESTINATION perf/bin)


