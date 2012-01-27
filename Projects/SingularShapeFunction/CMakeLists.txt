
########### next target ###############

file(GLOB header *.h)

set(main_SRCS main.cpp TExtFunction.cpp TDiscoFunction.cpp)

add_executable(SingularShapeFunction ${main_SRCS} ${header})

target_link_libraries(SingularShapeFunction pz pthread)

