# @file neopz/External/FAD/utils/CMakeLists.txt

target_include_directories(pz PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/External/FAD/utils>
        )

set(headers error.h promote.h timer.h tinyvec.h tinyveccpy.h vectors.h)
install(FILES ${headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/External/FAD/utils)
target_sources(pz PRIVATE ${headers})
