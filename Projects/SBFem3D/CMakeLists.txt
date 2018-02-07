file(GLOB headers *.h)
file(GLOB sources *.cpp)


set (MainRegular mainregular3dsolution.cpp Common3D.cpp ${headers})

add_executable(SBFem3DRegular3D ${MainRegular})

target_link_libraries(SBFem3DRegular3D pz)

