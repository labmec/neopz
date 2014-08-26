########### next target ###############

file(GLOB header *.h)

file(GLOB sources *.cpp)

add_executable(OilWaterSystem ${sources} ${header})

target_link_libraries(OilWaterSystem pz ${CBLAS_LIB} ${LAPACK_LIB}) 

file(COPY OilWaterLog4cxx.cfg DESTINATION .)
