
file(GLOB source *.cpp)
file(GLOB header *.h)

set(pztutormat ${pztutormat} ${source} ${header} PARENT_SCOPE)

add_executable(Perf02_SubStruct-Nathan-Timed ${source} ${header})

target_link_libraries(Perf02_SubStruct-Nathan-Timed pz pthread)

configure_file(8andares02.txt 8andares02.txt COPYONLY)



