
file(GLOB header *.h)

set(main_SRCS HydrFractPropag.cpp)

add_executable(HydraulicFracturePropagation ${main_SRCS} ${header})

target_link_libraries(HydraulicFracturePropagation pz pthread)

