# @file neopz/External/FAD/TinyFadET/CMakeLists.txt

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/External/FAD/TinyFadET>
                           )

set(headers tfad.h tfadfunc.h tfadlog.h tfadop.h)
install(FILES ${headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/External/FAD/TinyFadET)
target_sources(pz PRIVATE ${headers})
