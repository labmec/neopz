# @file neopz/Material/REAL/CMakeLists.txt  -- CMake file for REAL Materials

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
						   $<INSTALL_INTERFACE:.> 
						  )
file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE}/Material/REAL)

file(GLOB sources *.cpp)

if (USING_FAD AND BUILD_PLASTICITY_MATERIALS)
	if ((STATE_TYPE STREQUAL "double") AND (REAL_TYPE STREQUAL "double"))
		add_subdirectory(Plasticity)
	else()
		message ("Plasticity materials will not be built (set REAL_TYPE=double and STATE_TYPE=double to build them) !")
	endif()
endif()

set(pzmaterial_REAL_SRCS ${headers} ${sources} ${pzmaterial_Plasticity_SRCS} PARENT_SCOPE)
