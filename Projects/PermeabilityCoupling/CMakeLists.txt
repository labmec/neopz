
file(GLOB header *.h)

file(GLOB source *.cpp)

add_executable(PermeabilityCoupling ${source} ${header})

target_link_libraries(PermeabilityCoupling pz ${CBLAS_LIB} ${LAPACK_LIB})

