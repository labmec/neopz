
file(GLOB source *.cpp)
file(GLOB header *.h)

set(pztutormat ${pztutormat} ${source} ${header} PARENT_SCOPE)

add_executable(testmat_Tutorial ${source} ${header})

target_link_libraries(testmat_Tutorial pz pthread)

