# @file neopz/External/FAD/CMakeLists.txt  -- CMake file for the FAD external module

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/External/FAD>
                           )

set(headers fadType.h)

install(FILES ${headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/External/FAD)

target_sources(pz PRIVATE ${headers})

add_subdirectory(utils)
add_subdirectory(Fad)
add_subdirectory(TinyFad)
add_subdirectory(TinyFadET)
