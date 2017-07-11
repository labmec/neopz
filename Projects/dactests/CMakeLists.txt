
file(GLOB header *.h)

file(GLOB sourceEsfera mainEsfera.cpp )
add_executable(dacEsfera ${sourceEsfera} ${header})
target_link_libraries(dacEsfera pz)

file(GLOB sourceMeta PZMatPoissonD3.cpp LaplaceInCube.cpp LaplaceInQuadrilateral.cpp LaplaceInCircle.cpp LaplaceInCylinder.cpp LaplaceInSphere.cpp LaplaceInSolidSphere.cpp tools.cpp main.cpp)
add_executable(dacMeta ${sourceMeta} ${header})
target_link_libraries(dacMeta pz)

file(GLOB source2dPolinome main2dpolinomio.cpp PZMatPoissonD3.cpp TestHDivMesh.cpp )
add_executable(dacPolinome2D ${source2dPolinome} ${header})
target_link_libraries(dacPolinome2D pz)

file(GLOB source2dsin main2dseno.cpp PZMatPoissonD3.cpp TestHDivMesh.cpp )
add_executable(dacSin2D ${source2dsin} ${header})
target_link_libraries(dacSin2D pz)

file(GLOB source3dsenozc main3dzconst.cpp PZMatPoissonD3.cpp TestHDivMesh.cpp )
add_executable(dacZConst ${source3dsenozc} ${header})
target_link_libraries(dacZConst pz)

file(GLOB source3dsenoyc main3dyconst.cpp PZMatPoissonD3.cpp TestHDivMesh.cpp )
add_executable(dacYConst ${source3dsenoyc} ${header})
target_link_libraries(dacYConst pz)

file(GLOB source3dsenoxc main3dxconst.cpp PZMatPoissonD3.cpp TestHDivMesh.cpp )
add_executable(dacXConst ${source3dsenoxc} ${header})
target_link_libraries(dacXConst pz)

file(GLOB source3darctan main3darctan.cpp PZMatPoissonD3.cpp TestHDivMesh.cpp )
add_executable(dacArctan ${source3darctan} ${header})
target_link_libraries(dacArctan pz)

file(GLOB sourceAntigo PZMatPoissonD3.cpp TestHDivMesh.cpp main_antigo.cpp)
add_executable(dacAntigo ${sourceAntigo} ${header})
target_link_libraries(dacAntigo pz)