
if(BUILD_TUTORIAL)

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

	set(pzcurso ${pzcurso} ${pzclass1} ${pzclass2} ${pzclass3} ${pzclass4} ${pzclass5} PARENT_SCOPE )

endif()

