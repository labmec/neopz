file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(main_SRCS ${sources} ${headers})

#Plasticity works only with double type.
if ( STATE_TYPE STREQUAL "double" AND
     REAL_TYPE STREQUAL "double" AND 
     BUILD_PLASTICITY_MATERIALS )


	add_executable(PlasticityPVMC ${main_SRCS})
	target_link_libraries(PlasticityPVMC pz)

endif()