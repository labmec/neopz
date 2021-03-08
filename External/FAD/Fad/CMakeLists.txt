# @file neopz/External/FAD/Fad/CMakeLists.txt

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:${NEOPZ_INCLUDE}/External/FAD/Fad>
                          )

set(headers fad.h fadfunc.h fadlog.h fadop.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE}/External/FAD/Fad)
target_sources(pz PRIVATE ${headers})

