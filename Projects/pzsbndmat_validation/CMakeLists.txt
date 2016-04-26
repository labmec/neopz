########### next target ###############
file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(main_SRCS ${sources} ${headers})
add_executable(pzsbndmat_validation ${main_SRCS})
target_link_libraries(pzsbndmat_validation pz)