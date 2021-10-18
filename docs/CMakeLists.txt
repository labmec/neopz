#it will only come here if(BUILD_DOCS)

set(DOXYGEN_ADD_ALL OFF CACHE BOOL "Whether to add ALL source files to doxygen")

list(APPEND PZ_PUBLIC_HEADERS
    "${CMAKE_CURRENT_SOURCE_DIR}/pzmainpage.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/analysis/doxanalysis.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/structmatrix/doxstrmatrix.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/material/doxmaterial.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/solver/doxsolver.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/shape/doxshape.h"
    )

if(DOXYGEN_ADD_ALL)
    list(APPEND PZ_PUBLIC_HEADERS ${pz_sources})#this is valid just for this scope
endif()
string(REPLACE ";" " " DOXYGEN_INPUT_FILES "${PZ_PUBLIC_HEADERS}")
include(../cmake/EnableDoxygen.cmake)
enable_doxygen(pz)
if(BUILD_SPHINX_DOCS)
    if(NOT BUILD_DOCS)
        message(FATAL_ERROR "For building sphinx docu BUILD_DOCS needs to be enabled")
    endif()
    include(../cmake/EnableSphinx.cmake)
    enable_sphinx(pz)
endif()

set(docs_sources
    "${CMAKE_CURRENT_SOURCE_DIR}/pzmainpage.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in"
    "${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in"
    "${CMAKE_CURRENT_SOURCE_DIR}/index.rst"
    "${CMAKE_CURRENT_SOURCE_DIR}/references.bib")
add_subdirectory(material)
add_subdirectory(structmatrix)
add_subdirectory(solver)
    