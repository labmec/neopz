# @file neopz/External/FAD/TinyFad/CMakeLists.txt

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:${NEOPZ_INCLUDE}/External/FAD/TinyFad>
                          )

set(headers tinyfad.h tinyfadbin.h tinyfadfunc.h tinyfadlog.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE}/External/FAD/TinyFad)
target_sources(pz PRIVATE ${headers})

add_subdirectory(Specializations)
