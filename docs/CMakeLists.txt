#it will only come here if(BUILD_DOCS)

list(APPEND PZ_PUBLIC_HEADERS
    "${CMAKE_CURRENT_SOURCE_DIR}/pzmainpage.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/structmatrix/doxstrmatrix.h"
    )
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