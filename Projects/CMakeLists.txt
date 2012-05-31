# @file neopz/Projects/CMakeList.txt  -- First CMake file for Projects directory 

ENABLE_TESTING()
#
# COMPLEX NUMBER PROJECTS
#
# If your project use Complex numbers (STATE type = complex any type), please add it on the list below
if(BUILD_COMPLEX_PROJECTS)

	#Example to work with scalar as complex
	add_subdirectory(Helmholtz1D)

	#Examples used as tutorial for FEEC group (Prof Hugo)
	add_subdirectory(PDElinear1D)
	add_subdirectory(PDElinear2D)

else()

	if(BUILD_TUTORIAL)

		#Example to reading and using geometric mesh
		add_subdirectory(testgeom)

		#Testing and using shape functions
		add_subdirectory(testshape)

		#Curso PZ:
		add_subdirectory(CursoPZ)
		#add_subdirectory(PDElinear1D)
		#add_subdirectory(PDElinear2D)

		#Generic Example 2D with elasticity, using matrices phi, phi, etc
		add_subdirectory(testmat)

		#Using visual matrix of the DX
		add_subdirectory(visualmatrix)
		
		#Using discontinuous approximation functions
		add_subdirectory(descontinuo)

		#Doing h-refinement
	#	add_subdirectory(adaptivity)

		#source_group(CursoPZ FILES ${pzcurso})
		#source_group(Shape FILES ${pztutorshape})
		#source_group(Geometry FILES ${pztutorgeom})
		#source_group(Matrix FILES ${pztutormat})
		#source_group(VisualMatrixDX FILES ${visualmatrix})

		#set(pztutorial ${pztutorial} ${pzcurso} ${pztutorshape} ${pztutorgeom} ${pztutormat} ${visualmatrix})

		#add_custom_target(PZTutorial SOURCES ${pztutorial})


	endif(BUILD_TUTORIAL)

	if(BUILD_PROJECTS)
		
		add_subdirectory(MultipleLoad)
		# the auto adaptive process is working again!!!
		add_subdirectory(adaptivity)
		
		#Inserting project to hp adaptive development
		add_subdirectory(Adapt_HP_Jorge)

		add_subdirectory(Cubo_do_Nathan)
		
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
		add_subdirectory(QuadraticElements)
		add_subdirectory(ThermalStress)
		add_subdirectory(HDiv)
		add_subdirectory(AcoplamentoH1Hdiv)
		add_subdirectory(Validacao)
		add_subdirectory(Conductivity)
		add_subdirectory(ConsLaw)
		add_subdirectory(Descontinuo2)
		add_subdirectory(HydraulicFracturePropagation)

		add_subdirectory(SteamInjection)
		
		#projeto agnaldo 09-09-2011
		add_subdirectory(DoubleProjection1D)
		add_subdirectory(Poroelastic2D)
		add_subdirectory(SubStruct)
		add_subdirectory(Substruct_Viscoelastic)
		add_subdirectory(MultiphysicsSimulation)
		add_subdirectory(ProblemGirkmann)
		
		add_subdirectory(multplaca)

		#projeto Diogo 08-FEB-2012
		#add_subdirectory(Plasticity)
		#add_subdirectory(prodindex)

		#Projeto Teo Refinamento Paralelo
		add_subdirectory(RefThread)

		#Projeto Lua 2012
		add_subdirectory(Lua01)

		#Projeto Denise
		add_subdirectory(HdivPressureTeste)



		
		
	endif(BUILD_PROJECTS)


endif(BUILD_COMPLEX_PROJECTS)
