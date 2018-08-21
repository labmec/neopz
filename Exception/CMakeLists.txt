# @file neopz/Exception/CMakeLists.txt  -- CMake file for the Exception module

target_include_directories(pz PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:.> )

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE}/Exception)

file(GLOB sources *.cpp)

set(pzexception_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)

