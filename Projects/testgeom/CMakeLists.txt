########### next target ###############

set(main_SRCS geom.cpp)

add_executable(TestGeom ${main_SRCS})

configure_file(quad_st_800_1R_16X.gri quad_st_800_1R_16X.gri COPYONLY)

target_link_libraries(TestGeom pz pthread)

