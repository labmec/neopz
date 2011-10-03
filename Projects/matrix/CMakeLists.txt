########### next target ###############

file(GLOB header *.h)

set(main_SRCS main.cpp)

add_executable(matrix ${main_SRCS} ${header})

target_link_libraries(matrix pz pthread)

