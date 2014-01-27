include_directories(PlaneFracture)
include_directories(SWXGraph)

add_subdirectory(PlaneFracture)
add_subdirectory(SWXGraph)

file(GLOB header *.h)

set(main_SRCS HydrFractPropag.cpp)

add_executable(HydraulicFracturePropagation ${main_SRCS} ${header} ${pzplanefracture_STAT_SRCS} ${swxgraph_STAT_SRCS})

target_link_libraries(HydraulicFracturePropagation pz)

