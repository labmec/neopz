# @file neopz/External/FAD/TinyFad/CMakeLists.txt

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/External/FAD/TinyFad>
                           )

set(headers tinyfad.h tinyfadbin.h tinyfadfunc.h tinyfadlog.h)
install(FILES ${headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/External/FAD/TinyFad)
target_sources(pz PRIVATE ${headers})

add_subdirectory(Specializations)
