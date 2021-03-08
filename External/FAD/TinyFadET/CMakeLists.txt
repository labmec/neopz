# @file neopz/External/FAD/TinyFadET/CMakeLists.txt

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:${NEOPZ_INCLUDE}/External/FAD/TinyFadET>
                          )

set(headers tfad.h tfadfunc.h tfadlog.h tfadop.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE}/External/FAD/TinyFadET)
target_sources(pz PRIVATE ${headers})
