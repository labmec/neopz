file(GLOB source *.cpp)
file(GLOB header *.h)

set(pztutormat ${pztutormat} ${source} ${header} PARENT_SCOPE)

add_executable(SubStruct ${source} ${header})

target_link_libraries(SubStruct pz pthread)

configure_file(8andares02.txt 8andares02.txt COPYONLY)

configure_file(cube1.txt cube1.txt COPYONLY)
