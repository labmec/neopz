# @file neopz/Exception/CMakeLists.txt  -- CMake file for the Exception module

target_include_directories(pz PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/Exception>
                           )

set(headers
    TPZConvergenceException.h
    TPZException.h
    TPZInconsistentStateException.h
    )

set(sources
    TPZConvergenceException.cpp
    TPZException.cpp
    TPZInconsistentStateException.cpp
    )

install(FILES ${headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Exception)

target_sources(pz PRIVATE ${headers} ${sources})

