
if(BUILD_TUTORIAL)

	#Example using integration rules: IntegralExample
	add_subdirectory(cursoPZ_1)
	#Example using matrizes and solvers: MatrixExample
	add_subdirectory(cursoPZ_2)
	#Example using large mesh and topological elements (prismatic): TopologyExample
	add_subdirectory(cursoPZ_3)
	#Example testing refinement: RefineExample
	add_subdirectory(cursoPZ_4)
	#Example using writing and reading meshes: IOMeshExample
	add_subdirectory(cursoPZ_5)

endif()

