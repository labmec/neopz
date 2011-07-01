
if(BUILD_TUTORIAL)

	#add_subdirectory(Tutorial_Examples)
	#Testing and using shape functions
	add_subdirectory(testshape)
	#Example to reading and using geometric mesh
	add_subdirectory(testgeom)
	#Curso PZ:
		#Example using integration rules: IntegralExample
		#Example using matrizes and solvers: MatrixExample
		#Example using large mesh and topological elements (prismatic): TopologyExample
		#Example testing refinement: RefineExample
		#Example using writing and reading meshes: IOMeshExample
	add_subdirectory(cursoPZ)
		#add_subdirectory(cursoPZ_1)
		#add_subdirectory(cursoPZ_2)
		#add_subdirectory(cursoPZ_3)
		#add_subdirectory(cursoPZ_4)
		#add_subdirectory(cursoPZ_5)
	#Generic Example 2D with elasticity, using matrices phi, phi, etc
	add_subdirectory(testmat)
	#Using visual matrix of the DX
	add_subdirectory(visualmatrix)

endif()

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
	add_subdirectory(SteamInjection)
	add_subdirectory(ConsLaw)
	add_subdirectory(Descontinuo2)

if(BUILD_UNITTESTING)
	add_subdirectory(UnitTest)
endif()