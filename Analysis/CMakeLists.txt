# @file neopz/Analysis/CMakeLists.txt  -- CMake file for Analysis module. 

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:${NEOPZ_INCLUDE}/Analysis>
                          )

set(headers
    TPZNLMultGridAnalysis.h
    pzanalysis.h
    pzmganalysis.h
    pzsmanal.h
    pztransientanalysis.h
    doxanalysis.h
    pzanalysiserror.h
    pzeuleranalysis.h
    pznonlinanalysis.h
    pzsmfrontalanal.h
    )

set(sources
    TPZNLMultGridAnalysis.cpp
    pzanalysiserror.cpp
    pzeuleranalysis.cpp
    pznonlinanalysis.cpp
    pzsmfrontalanal.cpp
    tpzbiharmonicestimator.cpp
    pzanalysis.cpp
    pzblackoilanalysis.cpp
    pzmganalysis.cpp
    pzsmanal.cpp
    pztransientanalysis.cpp
    )

if (NOT BUILD_COMPLEX_PROJECTS)
    list(APPEND headers pzblackoilanalysis.h tpzbiharmonicestimator.h)
    list(APPEND sources pzblackoilanalysis.cpp tpzbiharmonicestimator.cpp)
endif ()

install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE}/Analysis)

target_sources(pz PRIVATE ${headers} ${sources})
