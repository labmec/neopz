
file(GLOB source *.cpp)
file(GLOB header *.h)

set(pztutorshape ${pztutorshape} ${source} ${header} PARENT_SCOPE)

add_executable(testshape_Tutorial ${source} ${header})

target_link_libraries(testshape_Tutorial pz pthread)

