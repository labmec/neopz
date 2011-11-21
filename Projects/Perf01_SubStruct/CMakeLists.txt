file(GLOB source *.cpp)
file(GLOB header *.h)

configure_file(8andares02.txt 8andares02.txt COPYONLY)
configure_file(log4cxx_perf.cfg log4cxx_perf.cfg COPYONLY)

set(pztutormat ${pztutormat} ${source} ${header} PARENT_SCOPE)

add_executable(Perf1_SubStruct ${source} ${header})

target_link_libraries(Perf1_SubStruct pz pthread)

