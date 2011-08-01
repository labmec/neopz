
file(GLOB source *.cpp)
file(GLOB header *.h)

set(pztutormat ${pztutormat} ${source} ${header} PARENT_SCOPE)

add_executable(TutorialMatrix ${source} ${header})

target_link_libraries(TutorialMatrix pz pthread)

