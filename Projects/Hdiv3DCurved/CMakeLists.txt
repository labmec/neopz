########### next target ###############

file(GLOB header *.h)

file(GLOB sources *.cpp)

add_executable(Hdiv3DCurved ${sources} ${header})

target_link_libraries(Hdiv3DCurved pz ${CBLAS_LIB} ${LAPACK_LIB}) 
