# @file neopz/External/Sloan/CMakeLists.txt  -- CMake file for the Sloan external module

set(headers sloan.h)
set(sources
        diamtr.cpp
        gegra.cpp
        iorele.cpp
        isorti.cpp
        label.cpp
        number.cpp
        profi1.cpp
        rootls.cpp
        vsrtp1.cpp
        )

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:${NEOPZ_INCLUDE}/External/Sloan>
                          )

install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE}/External/Sloan)

target_sources(pz PRIVATE ${headers} ${sources})