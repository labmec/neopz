########### next target ###############

file(GLOB header *.h)

file(GLOB sources *.cpp)

add_executable(ParallelPerformance ${sources} ${header})

target_link_libraries(ParallelPerformance pz)

