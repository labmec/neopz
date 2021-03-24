########### next target ###############

file(GLOB header *.h)

file(GLOB sources *.cpp)

add_executable(GradOfX ${sources} ${header})

target_link_libraries(GradOfX pz ${CBLAS_LIB} ${LAPACK_LIB}) 
