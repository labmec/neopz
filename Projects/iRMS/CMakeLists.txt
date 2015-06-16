########### next target ###############

file(GLOB header *.h)

file(GLOB sources *.cpp)

add_executable(iRMS ${sources} ${header})

target_link_libraries(iRMS pz ${CBLAS_LIB} ${LAPACK_LIB}) 

file(COPY iRMS.cfg DESTINATION .)

file(COPY RefPatterns.rpt DESTINATION .)
