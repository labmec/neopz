file(GLOB source *.cpp)
file(GLOB header *.h)

#set(pztutormat ${pztutormat} ${source} ${header} PARENT_SCOPE)

add_executable(Substruct_Viscoelastic ${source} ${header})

target_link_libraries(Substruct_Viscoelastic pz)

#configure_file(8andares02.txt 8andares02.txt COPYONLY)

#configure_file(cube1.txt cube1.txt COPYONLY)

include (CTestTestFile.cmake)