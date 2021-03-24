# @file neopz/External/CMakeLists.txt  -- CMake file for the External module

target_include_directories(pz PUBLIC
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/External>
                           )

set(public_headers
    TPZRenumbering.h
    )
set(headers
    TPZBoostGraph.h
    TPZCutHillMcKee.h
    TPZRenumbering.h
    TPZSloanRenumbering.h
    fblaswr.h
    pzsloan.h
    rcm.h
    tpznodesetcompute.h
    pzmetis.h
    )

set(sources
    TPZBoostGraph.cpp
    TPZCutHillMcKee.cpp
    TPZRenumbering.cpp
    TPZSloanRenumbering.cpp
    pzsloan.cpp
    rcm.cpp
    tpznodesetcompute.cpp
    pzmetis.cpp
    )

install(FILES ${public_headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/External)

target_sources(pz PRIVATE ${headers} ${sources})

add_subdirectory(sloan)

#for doxygen
foreach(header ${public_headers})
    list(APPEND PZ_PUBLIC_HEADERS " ${CMAKE_CURRENT_SOURCE_DIR}/${header}")
endforeach()
set(PZ_PUBLIC_HEADERS ${PZ_PUBLIC_HEADERS} PARENT_SCOPE)