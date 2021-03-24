# @file neopz/Integral/CMakeLists.txt  -- CMake file for the Integral module

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/Integral>
                           )

set(headers
    TPZIntQuadQuarterPoint.h
    adapt.h
    doxintegral.h
    pzquad.h
    tpzgaussrule.h
    tpzintpoints.h
    tpzintrulelist.h
    tpzintrulep3d.h
    tpzintrulet.h
    tpzintrulet3d.h
    tpzprinteg.h
    )

set(sources
    TPZIntQuadQuarterPoint.cpp
    adapt.cpp
    pzquad.cpp
    tpzgaussrule.cpp
    tpzintrulelist.cpp
    tpzintrulep3d.cpp
    tpzintrulet.cpp
    tpzintrulet3d.cpp
    tpzprinteg.cpp
    )

install(FILES ${headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Integral)

target_sources(pz PRIVATE ${headers} ${sources})
#for doxygen
foreach(header ${headers})
    list(APPEND PZ_PUBLIC_HEADERS " ${CMAKE_CURRENT_SOURCE_DIR}/${header}")
endforeach()
set(PZ_PUBLIC_HEADERS ${PZ_PUBLIC_HEADERS} PARENT_SCOPE)