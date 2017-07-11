
file(GLOB header *.h)

file(GLOB sources *.cpp)

add_executable(Aqp_EulerConsLaw ${sources} ${header} ${Euler_STAT_SRCS})

target_link_libraries(Aqp_EulerConsLaw pz)


