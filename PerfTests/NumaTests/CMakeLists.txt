file(GLOB source *.cpp)

file(GLOB header *.h)

add_executable(Perf-NumaTst ${source} ${header})
install(TARGETS Perf-NumaTst DESTINATION bin)
target_link_libraries(Perf-NumaTst pz Threads::Threads)

