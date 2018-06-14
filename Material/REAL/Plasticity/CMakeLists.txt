# @file neopz/Material/REAL/Plasticity/CMakeLists.txt  -- CMake file for Plasticity Materials

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
						   $<INSTALL_INTERFACE:.> 
						  )
file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE}/Material/REAL/Plasticity)

file(GLOB sources *.cpp)

set(pzmaterial_Plasticity_SRCS ${headers} ${sources} PARENT_SCOPE)
