# @file neopz/External/FAD/Fad/CMakeLists.txt

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/External/FAD/Fad>
                           )

set(headers fad.h fadfunc.h fadlog.h fadop.h)
install(FILES ${headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/External/FAD/Fad)
target_sources(pz PRIVATE ${headers})

