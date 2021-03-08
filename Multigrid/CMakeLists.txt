# @file neopz/Multigrid/CMakeLists.txt  -- CMake file for the multigrid module


target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:${NEOPZ_INCLUDE}/Multigrid>
                          )

set(headers pzmgsolver.h pztransfer.h)
set(sources pzmgsolver.cpp pztransfer.cpp)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE}/Multigrid)
target_sources(pz PRIVATE ${headers} ${sources})
