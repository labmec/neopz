
if(BUILD_TUTORIAL)
	#Testing and using shape functions
	add_subdirectory(testshape)

	#Example to reading and using geometric mesh
	add_subdirectory(testgeom)

	#Curso PZ:
	add_subdirectory(CursoPZ)

	#Generic Example 2D with elasticity, using matrices phi, phi, etc
	add_subdirectory(testmat)

	#Using visual matrix of the DX
	add_subdirectory(visualmatrix)
	
	source_group(CursoPZ FILES ${pzcurso})
	source_group(Shape FILES ${pztutorshape})
	source_group(Geometry FILES ${pztutorgeom})
	source_group(Matrix FILES ${pztutormat})
	source_group(VisualMatrixDX FILES ${visualmatrix})

	set(pztutorial ${pztutorial} ${pzcurso} ${pztutorshape} ${pztutorgeom} ${pztutormat} ${visualmatrix})

	add_custom_target(PZTutorial SOURCES ${pztutorial})


endif(BUILD_TUTORIAL)

if(BUILD_PROJECTS)
	add_subdirectory(HyperElasticityFAD)

	#Increment singular shape function: SingularShapeFunction
	add_subdirectory(SingularShapeFunction)
	add_subdirectory(angletest)
	add_subdirectory(checkconv)

	#Orthogonalization process of Gram-Schmidt
	add_subdirectory(matrix)

	#Inserting material and boundary conditions by function
	add_subdirectory(termica)
	add_subdirectory(Steam)
	add_subdirectory(ThermalStress)
	add_subdirectory(HDiv)
	add_subdirectory(Conductivity)
	add_subdirectory(ConsLaw)
	add_subdirectory(Descontinuo2)

	add_subdirectory(SteamInjection)

endif(BUILD_PROJECTS)
