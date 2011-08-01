
file(GLOB source *.cpp)
file(GLOB header *.h)

set(pztutorshape ${pztutorshape} ${source} ${header} PARENT_SCOPE)

add_executable(TutorialShape ${source} ${header})

target_link_libraries(TutorialShape pz pthread)

