# @file neopz/External/CMakeLists.txt  -- CMake file for the External module

target_include_directories(pz PUBLIC
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/External>
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
    )

set(sources
    TPZBoostGraph.cpp
    TPZCutHillMcKee.cpp
    TPZRenumbering.cpp
    TPZSloanRenumbering.cpp
    pzsloan.cpp
    rcm.cpp
    tpznodesetcompute.cpp
    )

if(USING_METIS)
    LIST(APPEND headers pzmetis.h)
    LIST(APPEND sources pzmetis.cpp)
endif()

install(FILES ${headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/External)

target_sources(pz PRIVATE ${headers} ${sources})

add_subdirectory(sloan)

#enabling FAD library
add_subdirectory(FAD)
