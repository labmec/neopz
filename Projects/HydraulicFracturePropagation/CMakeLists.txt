
set(main_SRCS HydrFractPropag.cpp)

add_executable(HydraulicFracturePropagation ${main_SRCS})

target_link_libraries(HydraulicFracturePropagation pz pthread)

