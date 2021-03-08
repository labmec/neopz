file(GLOB source *.cpp)

file(GLOB header *.h)

add_executable(migration-perf ${source} ${header})
install(TARGETS migration-perf DESTINATION bin)
target_link_libraries(migration-perf pz Threads::Threads)

