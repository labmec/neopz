# @file neopz/External/FAD/CMakeLists.txt  -- CMake file for the FAD external module

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
						   $<INSTALL_INTERFACE:.> 
						  )

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE}/External/FAD)

file(GLOB utils_headers utils/*.h)
install(FILES ${utils_headers} DESTINATION ${NEOPZ_INCLUDE}/External/FAD/utils)

file(GLOB sources *.cpp)

add_subdirectory(Fad)
add_subdirectory(TinyFadET)
add_subdirectory(TinyFad)

set(FAD_fad_SRCS ${FAD_fad_SRCS} PARENT_SCOPE)

set(FAD_Tiny_SRCS ${FAD_Tiny_SRCS} PARENT_SCOPE)

set(FAD_TinyET_SRCS ${FAD_TinyET_SRCS} PARENT_SCOPE)

set(FAD_SRCS ${headers} ${sources} PARENT_SCOPE)

