
target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:${NEOPZ_INCLUDE}/Optimization>
                          )

set(headers TPZStochasticSearch.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE}/Optimization)
target_sources(pz PRIVATE ${header})
