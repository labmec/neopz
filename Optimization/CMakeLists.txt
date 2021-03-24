
target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/Optimization>
                           )

set(headers TPZStochasticSearch.h)
install(FILES ${headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Optimization)
target_sources(pz PRIVATE ${header})

#for doxygen
foreach(header ${headers})
    list(APPEND PZ_PUBLIC_HEADERS " ${CMAKE_CURRENT_SOURCE_DIR}/${header}")
endforeach()
set(PZ_PUBLIC_HEADERS ${PZ_PUBLIC_HEADERS} PARENT_SCOPE)