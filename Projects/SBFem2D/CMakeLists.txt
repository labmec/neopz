file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(MymainOLD Main.cpp MainJSon.cpp Common.cpp ${headers})

add_executable(SBFem2D ${MymainOLD})

target_link_libraries(SBFem2D pz)

set (MainRegular mainregularsolution.cpp Common.cpp ${headers})

add_executable(SBFem2DRegular ${MainRegular})

target_link_libraries(SBFem2DRegular pz)

set (MainParabolic main2DParabolicEq.cpp Common.cpp ${headers})

add_executable(SBFem2DParabolic ${MainParabolic})

target_link_libraries(SBFem2DParabolic pz)

set (MainCircle maincircle.cpp Common.cpp ${headers})

add_executable(SBFem2DCircle ${MainCircle})

target_link_libraries(SBFem2DCircle pz)

set (MainHeterogeneous mainheterogeneous.cpp Common.cpp ${headers})

add_executable(SBFem2DHeterogeneous ${MainHeterogeneous})

target_link_libraries(SBFem2DHeterogeneous pz)

set (MainCheckerboard maincheckerboard.cpp Common.cpp MainJSon.cpp ${headers})

add_executable(SBFem2DCheckerboard ${MainCheckerboard})

target_link_libraries(SBFem2DCheckerboard pz)

# SBFem2D with scaling center at the boundary of the domain and Dirichlet imposed
set (MainSqrt mainsqrt.cpp Common.cpp ${headers})

add_executable(SBFem2DSqrt ${MainSqrt})

target_link_libraries(SBFem2DSqrt pz)

