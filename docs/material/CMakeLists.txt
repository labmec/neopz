list(APPEND docs_sources
    "${CMAKE_CURRENT_SOURCE_DIR}/doxmaterial.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/index.rst"
    "${CMAKE_CURRENT_SOURCE_DIR}/availablemats.rst"
    "${CMAKE_CURRENT_SOURCE_DIR}/combinedspaces.rst"
    "${CMAKE_CURRENT_SOURCE_DIR}/elasticity.rst "
    "${CMAKE_CURRENT_SOURCE_DIR}/electromagnetics.rst"
    "${CMAKE_CURRENT_SOURCE_DIR}/lagrangemat.rst"
    "${CMAKE_CURRENT_SOURCE_DIR}/nullmat.rst"
    "${CMAKE_CURRENT_SOURCE_DIR}/poisson.rst"
    "${CMAKE_CURRENT_SOURCE_DIR}/projections.rst"
    "${CMAKE_CURRENT_SOURCE_DIR}/singlespace.rst"
    )

set(docs_sources ${docs_sources} PARENT_SCOPE)