
file(GLOB header *.h)

file(GLOB source *.cpp)

add_executable(dactests ${source} ${header})

target_link_libraries(dactests pz ${CBLAS_LIB} ${LAPACK_LIB})

