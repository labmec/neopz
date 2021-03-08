# @file neopz/Util/Hash/CMakeLists.txt  -- CMake file for the Hash module

target_include_directories(pz PUBLIC
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:${NEOPZ_INCLUDE}/Util/Hash>
                           )

set(headers MurmurHash3.h TPZHash.h)
set(sources MurmurHash3.cpp TPZHash.cpp)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE}/Util/Hash)
target_sources(pz PRIVATE ${headers} ${sources})
