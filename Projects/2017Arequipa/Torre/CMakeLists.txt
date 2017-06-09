file(GLOB header *.h)

file(GLOB sources *.cpp)

add_executable(Aqp_SWXTorre ${sources} ${header} ${swxLib_STAT_SRCS})

target_link_libraries(Aqp_SWXTorre pz)


