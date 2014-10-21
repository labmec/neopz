
file(GLOB header *.h)

file(GLOB sourceMeta PZMatPoissonD3.cpp TestHDivMesh.cpp main.cpp)
add_executable(dacMeta ${sourceMeta} ${header})
target_link_libraries(dacMeta pz ${CBLAS_LIB} ${LAPACK_LIB})

file(GLOB source2dPolinome main2dpolinomio.cpp PZMatPoissonD3.cpp TestHDivMesh.cpp )
add_executable(dacPolinome2D ${source2dPolinome} ${header})
target_link_libraries(dacPolinome2D pz ${CBLAS_LIB} ${LAPACK_LIB})

file(GLOB source2dPolinome main3dzconst.cpp PZMatPoissonD3.cpp TestHDivMesh.cpp )
add_executable(dacZConst ${source2dPolinome} ${header})
target_link_libraries(dacZConst pz ${CBLAS_LIB} ${LAPACK_LIB})



