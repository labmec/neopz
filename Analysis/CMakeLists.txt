# @file neopz/Analysis/CMakeLists.txt  -- CMake file for Analysis module. 

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/Analysis>
                           )
set(public_headers
  pzanalysis.h
  pznonlinanalysis.h
  pzsmanal.h #TPZSubMeshAnalysis
  )
set(headers
    TPZNLMultGridAnalysis.h
    pzanalysis.h
    pzmganalysis.h
    pzsmanal.h
    pztransientanalysis.h
    doxanalysis.h
    pzeuleranalysis.h
    pznonlinanalysis.h
    pzsmfrontalanal.h
    )

set(sources
    TPZNLMultGridAnalysis.cpp
    pzeuleranalysis.cpp
    pznonlinanalysis.cpp
    pzsmfrontalanal.cpp
    pzanalysis.cpp
    pzblackoilanalysis.cpp
    pzmganalysis.cpp
    pzsmanal.cpp
    pztransientanalysis.cpp
    )

if (NOT BUILD_COMPLEX_PROJECTS)
    list(APPEND headers pzblackoilanalysis.h)
    list(APPEND sources pzblackoilanalysis.cpp)
endif ()

install(FILES ${public_headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Analysis)

target_sources(pz PRIVATE ${headers} ${sources})

#for doxygen
foreach(header ${public_headers})
    list(APPEND PZ_PUBLIC_HEADERS " ${CMAKE_CURRENT_SOURCE_DIR}/${header}")
endforeach()
set(PZ_PUBLIC_HEADERS ${PZ_PUBLIC_HEADERS} PARENT_SCOPE)