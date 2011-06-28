
########### next target ###############

set(main_SRCS main.cpp TExtFunction.cpp TDiscoFunction.cpp)

add_executable(discShape ${main_SRCS})

target_link_libraries(discShape pz pthread)

