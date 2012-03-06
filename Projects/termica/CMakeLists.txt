########### next target ###############

file(GLOB header *.h)

set(main_SRCS termi.cpp fast.cpp)

configure_file(quad_st_800_1R_16X.gri quad_st_800_1R_16X.gri COPYONLY)
#configure_file(geomesh_Crude.txt geomesh_Crude.txt COPYONLY)
#configure_file(geomesh_Fine.txt geomesh_Fine.txt COPYONLY)


add_executable(termica ${main_SRCS} ${header})

target_link_libraries(termica pz pthread)

