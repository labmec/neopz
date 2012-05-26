file(GLOB source *.cpp)

file(GLOB header *.h)

add_executable(Perf-SubStruct ${source} ${header})

target_link_libraries(Perf-SubStruct pz pthread)

