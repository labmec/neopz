file(GLOB headers *.h ../HydraulicFracturePropagation/PlaneFracture/TPZJIntegral.h)
file(GLOB sources *.cpp ../HydraulicFracturePropagation/PlaneFracture/TPZJIntegral.cpp)

set(main_SRCS ${sources} ${headers})

add_executable(ReducedSpace ${main_SRCS} ${headers})

target_link_libraries(ReducedSpace pz)


