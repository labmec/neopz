########### next target ###############

file(GLOB header *.h)

file(GLOB sources *.cpp)

add_executable(ThermalStress ${sources} ${header})

target_link_libraries(ThermalStress pz pthread)

