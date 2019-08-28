
file(GLOB header *.h)

file(GLOB source *.cpp)

add_executable(GeoMechanicROM ${source} ${header})

target_link_libraries(GeoMechanicROM pz ${CBLAS_LIB} ${LAPACK_LIB})

