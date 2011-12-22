
file(GLOB source *.cpp)
file(GLOB header *.h)

set(pztutormat ${pztutormat} ${source} ${header} PARENT_SCOPE)

add_executable(Perf02_SubStruct-Merge ${source} ${header})

target_link_libraries(Perf02_SubStruct-Merge pz pthread)

configure_file(8andares02.txt 8andares02.txt COPYONLY)



