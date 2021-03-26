set(main substruct.cpp)
set(test substruct_test.cpp)

#set(pztutormat ${pztutormat} ${source} ${header} PARENT_SCOPE)

add_executable(SubStruct ${main})

add_executable(SubStruct_Test ${test})

target_link_libraries(SubStruct pz)

target_link_libraries(SubStruct_Test pz)

#configure_file(8andares02.txt 8andares02.txt COPYONLY)

#configure_file(cube1.txt cube1.txt COPYONLY)

if(USING_METIS)
include (CTestTestFile.cmake)
endif()