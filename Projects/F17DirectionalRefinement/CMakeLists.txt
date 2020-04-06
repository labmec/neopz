

########### next target ###############

set(F17_SRCS
    F17.cpp)

add_executable(F17DirectionalRef ${F17_SRCS})

target_link_libraries(F17DirectionalRef pz)
