########### next target ###############

file(GLOB header *.h)

file(GLOB sources *.cpp)

add_executable(OilWaterSystem ${sources} ${header})

target_link_libraries(OilWaterSystem pz)

