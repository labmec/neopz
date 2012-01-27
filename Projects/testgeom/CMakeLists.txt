# @file neopz/Projects/testgeom/CMakeLists.txt  -- CMake file for the geometry tutorial  

########### next target ###############

file(GLOB source *.cpp)
file(GLOB header *.h)

set(pztutorgeom ${pztutorgeom} ${source} ${header} PARENT_SCOPE)

add_executable(testgeom_TutorialGeo geom.cpp pzgbmesh.cpp ${header})

configure_file(quad_st_800_1R_16X.gri quad_st_800_1R_16X.gri COPYONLY)

target_link_libraries(testgeom_TutorialGeo pz pthread)


add_executable(testgeom_TutorialComp pzcompmesh.cpp ${header})

target_link_libraries(testgeom_TutorialComp pz pthread)
