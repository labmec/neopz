file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(main_SRCS ${sources} ${headers})

if(STATE_COMPLEX)

	add_executable(Helmholtz1D ${main_SRCS})

	target_link_libraries(Helmholtz1D pz )

else()

	add_executable(HelmholtzComplex1D ${main_SRCS})

	target_link_libraries(HelmholtzComplex1D pz)

endif()