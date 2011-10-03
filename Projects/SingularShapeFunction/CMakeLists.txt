
########### next target ###############

file(GLOB header *.h)

set(main_SRCS main.cpp TExtFunction.cpp TDiscoFunction.cpp)

add_executable(discShape ${main_SRCS} ${header})

target_link_libraries(discShape pz pthread)

