# @file neopz/Material/Complex/CMakeLists.txt  -- CMake file for Complex Materials


target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
			   $<INSTALL_INTERFACE:${NEOPZ_INCLUDE}/Material/Complex> 
						  )

set(headers "")
set(sources "")

install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE}/Material/Complex)
target_sources(pz PRIVATE ${headers} ${sources})