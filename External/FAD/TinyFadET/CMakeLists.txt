# @file neopz/External/FAD/TinyFadET/CMakeLists.txt

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
						   $<INSTALL_INTERFACE:.> 
						  )

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE}/External/FAD/TinyFadET/)

file(GLOB sources *.cpp)

set(FAD_TinyET_SRCS ${headers} ${sources} PARENT_SCOPE)
