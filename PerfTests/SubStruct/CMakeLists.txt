file(GLOB source *.cpp)

file(GLOB header *.h)

add_executable(Perf-SubStruct ${source} ${header})
install(TARGETS Perf-SubStruct DESTINATION bin)
target_link_libraries(Perf-SubStruct pz)

configure_file(8andares02.txt 8andares02.txt COPYONLY)
configure_file(cube1.txt cube1.txt COPYONLY)

