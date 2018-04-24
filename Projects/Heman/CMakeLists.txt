

########### next target ###############

set(HeMan_SRCS
    c0-simplequad.cpp
    c1-simplecube.cpp
    F17.cc
    MeshReader.cpp
    Quad_1_Rib_Marked.cpp
    Quad_2_Rib_Ad_Marked.cpp)

add_executable(Heman ${HeMan_SRCS})

target_link_libraries(Heman pz)

set(RefPoco_SRCS
    RefPoco.cpp)

add_executable(GeracaoPoco ${RefPoco_SRCS})

target_link_libraries(GeracaoPoco pz)
