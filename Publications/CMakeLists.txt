# @file neopz/Publications/CMakeLists.txt  -- CMake file for the Publications module

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
						   $<INSTALL_INTERFACE:.> 
						  )

if(NOT BUILD_COMPLEX_PROJECTS) ##projetos nao compativeis com state=complex<double>
file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE}/Publications)

file(GLOB sources *.cpp)

set(pzpublications_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)
endif()
