
if(BUILD_TUTORIAL)

	#Example using gradient reconstruction by element based on least squares
	add_subdirectory(GradientReconstructed)

	#Example using integration rules: IntegralExample
	add_subdirectory(Class_1)
	#Example using matrizes and solvers: MatrixExample
	add_subdirectory(Class_2)
	#Example using large mesh and topological elements (prismatic): TopologyExample
	add_subdirectory(Class_3)
	#Example testing refinement: RefineExample
	add_subdirectory(Class_4)
	#Example using writing and reading meshes: IOMeshExample
	add_subdirectory(Class_5)
	#Example : Elasticity Material 2D
	add_subdirectory(Class_6)

	#Example : Elasticity Material 2D
	add_subdirectory(MiProyecto)
	#Example : Elasticity Material 2D
	add_subdirectory(Problema1D)

	#Example : 
	add_subdirectory(ExemploSimples)
	#Example : 
	add_subdirectory(ExemploSimples2D)
	#Example : Geometric Mapping
	add_subdirectory(ExemploMapeamentoGeometrico)


	#Example made to LNCC course (Simworx - Pira and Labmec - Phil)
	add_subdirectory(Problema2D)

	set(pzcurso ${pzcurso} ${pzclass1} ${pzclass2} ${pzclass3} ${pzclass4} ${pzclass5} ${pzclass6} PARENT_SCOPE )

endif()

