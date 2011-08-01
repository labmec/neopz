########### next target ###############

file(GLOB source *.cpp)
file(GLOB header *.h)

set(pztutorgeom ${pztutorgeom} ${source} ${header} PARENT_SCOPE)

add_executable(TutorialGeometry ${source} ${header})

configure_file(quad_st_800_1R_16X.gri quad_st_800_1R_16X.gri COPYONLY)

target_link_libraries(TutorialGeometry pz pthread)

