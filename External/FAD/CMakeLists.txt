# @file neopz/External/FAD/CMakeLists.txt  -- CMake file for the FAD external module

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:${NEOPZ_INCLUDE}/External/FAD>
                          )

set(headers fadType.h)

install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE}/External/FAD)

target_sources(pz PRIVATE ${headers})

add_subdirectory(utils)
add_subdirectory(Fad)
add_subdirectory(TinyFadET)
add_subdirectory(TinyFad)
