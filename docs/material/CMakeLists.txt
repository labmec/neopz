list(APPEND docs_sources
    "${CMAKE_CURRENT_SOURCE_DIR}/doxmaterial.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/index.rst"
    "${CMAKE_CURRENT_SOURCE_DIR}/availablemats.rst"
    "${CMAKE_CURRENT_SOURCE_DIR}/singlespace.rst"
    "${CMAKE_CURRENT_SOURCE_DIR}/commoninterfaces.rst"
    "${CMAKE_CURRENT_SOURCE_DIR}/combinedspaces.rst"
    "${CMAKE_CURRENT_SOURCE_DIR}/elasticity.inc "
    "${CMAKE_CURRENT_SOURCE_DIR}/electromagnetics.inc"
    "${CMAKE_CURRENT_SOURCE_DIR}/lagrangemat.inc"
    "${CMAKE_CURRENT_SOURCE_DIR}/nullmat.inc"
    "${CMAKE_CURRENT_SOURCE_DIR}/poisson.inc"
    "${CMAKE_CURRENT_SOURCE_DIR}/projections.inc"
    )

set(docs_sources ${docs_sources} PARENT_SCOPE)