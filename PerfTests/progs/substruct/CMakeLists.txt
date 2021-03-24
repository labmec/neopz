file(GLOB source *.cpp)
file(GLOB header *.h)

add_executable(substruct-perf ${source} ${header})
target_link_libraries(substruct-perf pz)
install(TARGETS substruct-perf DESTINATION perf/bin)

