file(GLOB source *.cpp)

file(GLOB header *.h)

add_executable(Perf-TBB ${source} ${header})
install(TARGETS Perf-TBB DESTINATION bin)
target_link_libraries(Perf-TBB pz)

