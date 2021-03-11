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
    pzmetis.h
    pzsloan.h
    rcm.h
    tpznodesetcompute.h
    )

set(sources
    TPZBoostGraph.cpp
    TPZCutHillMcKee.cpp
    TPZRenumbering.cpp
    TPZSloanRenumbering.cpp
    pzmetis.cpp
    pzsloan.cpp
    rcm.cpp
    tpznodesetcompute.cpp
    )

install(FILES ${headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/External)

target_sources(pz PRIVATE ${headers} ${sources})

add_subdirectory(sloan)

#enabling FAD library
option(USING_FAD "Whether the FAD package will be linked with the library" OFF)
if(USING_FAD)
	add_subdirectory(FAD)
endif(USING_FAD)
