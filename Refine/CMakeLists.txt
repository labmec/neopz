# @file neopz/Refine/CMakeLists.txt  -- CMake file for the Refine module

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/Refine>
                           )

set(headers
    TPZRefCube.h
    TPZRefLinear.h
    TPZRefPattern.h
    TPZRefPatternDataBase.h
    TPZRefPatternTools.h
    doxrefine.h
    pzrefpoint.h
    pzrefprism.h
    pzrefpyram.h
    pzrefquad.h
    pzreftetrahedra.h
    pzreftriangle.h
   )

set(sources
    TPZRefCube.cpp
    TPZRefPattern.cpp
    TPZRefPatternTools.cpp
    pzrefprism.cpp
    pzrefquad.cpp
    pzreftriangle.cpp
    TPZRefLinear.cpp
    TPZRefPatternDataBase.cpp
    pzrefpoint.cpp
    pzrefpyram.cpp
    pzreftetrahedra.cpp
   )

install(FILES ${headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Refine)

target_sources(pz PRIVATE ${headers} ${sources})

# files with description of the refinement patterns
add_subdirectory(RefPatterns)

#for doxygen
foreach(header ${headers})
    list(APPEND PZ_PUBLIC_HEADERS " ${CMAKE_CURRENT_SOURCE_DIR}/${header}")
endforeach()
set(PZ_PUBLIC_HEADERS ${PZ_PUBLIC_HEADERS} PARENT_SCOPE)