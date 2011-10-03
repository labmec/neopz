########### next target ###############

file(GLOB header *.h)

set(main_SRCS termi.cpp fast.cpp)

add_executable(termica ${main_SRCS} ${header})

target_link_libraries(termica pz pthread)

