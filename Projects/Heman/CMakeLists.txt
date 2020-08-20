

########### next target ###############

set(RefPoco_SRCS
    RefPoco.cpp)

add_executable(GeracaoPoco ${RefPoco_SRCS})

target_link_libraries(GeracaoPoco pz)
