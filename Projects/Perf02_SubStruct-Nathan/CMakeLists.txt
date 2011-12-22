
file(GLOB source *.cpp)
file(GLOB header *.h)

set(pztutormat ${pztutormat} ${source} ${header} PARENT_SCOPE)

add_executable(Perf02_SubStruct-Nathan ${source} ${header})

target_link_libraries(Perf02_SubStruct-Nathan pz pthread)

configure_file(8andares02.txt 8andares02.txt COPYONLY)



