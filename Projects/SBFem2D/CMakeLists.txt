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

set (MainStressFiles mainstressintensity.cpp Common.cpp ${headers})

add_executable(SBFem2DStressIntensity ${MainStressFiles})

target_link_libraries(SBFem2DStressIntensity pz)

set (MainCrack2DFiles mainEmbeddedFractures.cpp Common.cpp Common3D.cpp ${headers})

add_executable(SBFem2DEmbeddedFractures ${MainCrack2DFiles})

target_link_libraries(SBFem2DEmbeddedFractures pz)

configure_file(TwoOrthogonalCracks.txt TwoOrthogonalCracks.txt COPYONLY)

set (MainSBFemFemFiles mainsbfem_fem.cpp Common.cpp ${headers})

add_executable(SBFem2DFem ${MainSBFemFemFiles})

target_link_libraries(SBFem2DFem pz)

