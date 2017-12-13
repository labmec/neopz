file(GLOB headers *.h)
file(GLOB sources *.cpp)


set (MainRegular mainregular3dsolution.cpp Common3D.cpp ${headers})

add_executable(SBFem3DRegular3D ${MainRegular})

target_link_libraries(SBFem3DRegular3D pz)

set (MainCook maincooksimulation.cpp Common3D.cpp ${headers})

add_executable(SBFem3DCook ${MainCook})

target_link_libraries(SBFem3DCook pz)

set (MainSphere mainspheresimulation.cpp Common3D.cpp ${headers})

add_executable(SBFem3DSphere ${MainSphere})

target_link_libraries(SBFem3DSphere pz)

configure_file(CooksMembrane_poly_16_1_1.txt CooksMembrane_poly_16_1_1.txt COPYONLY)

configure_file(CooksMembrane_sbfemesh_32_2_1.txt CooksMembrane_sbfemesh_32_2_1.txt COPYONLY)

configure_file(CooksMembrane_sbfemesh_64_2_1.txt CooksMembrane_sbfemesh_64_2_1.txt COPYONLY)

configure_file(CooksMembrane_sbfemesh_128_4_1.txt CooksMembrane_sbfemesh_128_4_1.txt COPYONLY)

configure_file(spheres_10_50_sbfemesh_64_8_1.txt spheres_10_50_sbfemesh_64_8_1.txt COPYONLY)

configure_file(spheres_10_50_sbfemesh_32_8_1.txt spheres_10_50_sbfemesh_32_8_1.txt COPYONLY)

configure_file(spheres_10_50_sbfemesh_128_8_1.txt spheres_10_50_sbfemesh_128_8_1.txt COPYONLY)