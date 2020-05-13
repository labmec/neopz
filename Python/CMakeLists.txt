set(CMAKE_INCLUDE_CURRENT_DIR ON)

include_directories(${PROJECT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}
                   "${PROJECT_SOURCE_DIR}/Analysis"
                   "${PROJECT_SOURCE_DIR}/Common"
                   "${PROJECT_BINARY_DIR}/Common"                    # pz_config.h
                   "${PROJECT_SOURCE_DIR}/External/"
                   "${PROJECT_SOURCE_DIR}/Frontal/"
                   "${PROJECT_SOURCE_DIR}/Geom"
                   "${PROJECT_SOURCE_DIR}/Integral"
                   "${PROJECT_SOURCE_DIR}/Material"
                   "${PROJECT_SOURCE_DIR}/Material/REAL"
                   "${PROJECT_SOURCE_DIR}/Matrix"
                   "${PROJECT_SOURCE_DIR}/Mesh"
                   "${PROJECT_SOURCE_DIR}/Multigrid"
                   "${PROJECT_SOURCE_DIR}/Pre"
                   "${PROJECT_SOURCE_DIR}/Post"
                   "${PROJECT_SOURCE_DIR}/Refine"
                   "${PROJECT_SOURCE_DIR}/Save"
                   "${PROJECT_SOURCE_DIR}/SpecialMaps"
                   "${PROJECT_SOURCE_DIR}/StrMatrix"
                   "${PROJECT_SOURCE_DIR}/Topology"
                   "${PROJECT_SOURCE_DIR}/Util"
                    ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_SOURCE_DIR} )

if(USING_FAD)
    include_directories( "${PROJECT_SOURCE_DIR}/External/FAD/"
                         "${PROJECT_SOURCE_DIR}/External/FAD/Fad"
                         "${PROJECT_SOURCE_DIR}/External/FAD/TinyFad"
                         "${PROJECT_SOURCE_DIR}/External/FAD/TinyFadET")
endif(USING_FAD)

if(BUILD_PLASTICITY_MATERIALS)
    include_directories( "${PROJECT_SOURCE_DIR}/Material/REAL/Plasticity" )
endif()


# pybind11
find_package(pybind11 REQUIRED)

# Project binding sources
pybind11_add_module( NEOPZ src/neopz_bindings.cpp )

IF(APPLE)
  target_link_libraries(NEOPZ PRIVATE pz)
ELSE()
  IF(${USING_MKL})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-narrowing -Wl,--no-as-needed")
    find_library(MKL_LIB_DEF NAMES libmkl_def.so libmkl_def.dylib PATHS ${SEARCH_DIRS} /opt/intel/mkl/lib /opt/intel/mkl/lib/intel64/ /softwares/intel/mkl/lib /softwares/intel/mkl/lib/intel64/)
    target_link_libraries(NEOPZ PUBLIC pz ${MKL_LIB_DEF})
    target_include_directories(NEOPZ PUBLIC ${MKL_INCLUDE})
  ELSE()
    target_link_libraries(NEOPZ PUBLIC pz)
  ENDIF()
ENDIF()

