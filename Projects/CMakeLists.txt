# @file neopz/Projects/CMakeList.txt  -- First CMake file for Projects directory

ENABLE_TESTING()

	if(BUILD_TUTORIAL)

		if(BUILD_COMPLEX_PROJECTS)
			#add_subdirectory(PDElinear1D)
			#add_subdirectory(PDElinear2D)
		else(BUILD_COMPLEX_PROJECTS)

			source_group(CursoPZ FILES ${pzcurso})
			#source_group(Shape FILES ${pztutorshape})
			#source_group(Geometry FILES ${pztutorgeom})
			#source_group(Matrix FILES ${pztutormat})
			#source_group(VisualMatrixDX FILES ${visualmatrix})

			#Curso em Arequipa 2017
			#source_group(Aqp_2017 FILES ${aqp2017})
			#add_subdirectory(2017Arequipa)

			#set(pztutorial ${pztutorial} ${pzcurso} ${pztutorshape} ${pztutorgeom} ${pztutormat} ${visualmatrix})

			#add_custom_target(PZTutorial SOURCES ${pztutorial})
			#add_subdirectory(Finite3DElasticity)

			#Curso PZ:
			add_subdirectory(CursoPZ)

			#Example to reading and using geometric mesh
			#add_subdirectory(testgeom)

			#Testing and using shape functions
			#add_subdirectory(testshape)

			#Generic Example 2D with elasticity and 2D linear problem showing an example call to Contribute
			#add_subdirectory(testmat)

			#Using visual matrix of the DX
			add_subdirectory(visualmatrix)


			# project to save and load a computational mesh
			add_subdirectory(SaveAndLoad)

			# Tutorial about hp-adaptive techniques
			add_subdirectory(Poisson3D_Fichera)
			add_subdirectory(Poisson3D_Shock)
			# the auto adaptive process is working again!!!
			#add_subdirectory(adaptivity)

			#Example of the  multi-physics simulation with discontinuous Galerkin
 			add_subdirectory(PauloBosing)

			#Example of the tracer flow problem in multiphysics simulation.
  			add_subdirectory(TracadorPassivo)

			#Project to test the method hdiv stabilized.
			add_subdirectory(HDivEstabilizado)
			add_subdirectory(HDivEstabHHetero)


			#Example of the hybrid mesh
			add_subdirectory(HybridMesh)

			#Example of the  convection problem
			add_subdirectory(Problema_Transporte)


			# Simple project to test rational function defined by Cedric on pyramids
			add_subdirectory(RationalFunctionPyramid)

			add_subdirectory(proDG)

		endif(BUILD_COMPLEX_PROJECTS)


	endif(BUILD_TUTORIAL)

	if(BUILD_BROKEN)
		if(BUILD_COMPLEX_PROJECTS)
			#add_subdirectory(PDElinear1D)
			#add_subdirectory(PDElinear2D)
		else(BUILD_COMPLEX_PROJECTS)
                	#Testing gradient reconstruction
			add_subdirectory(Gradient_Reconstruction)
			#exercise on gradient reconstruction, but the results are wrong!
			add_subdirectory(Poisson2D_ArcTan_Disc)


		endif(BUILD_COMPLEX_PROJECTS)
	endif(BUILD_BROKEN)

	if(BUILD_PROJECTS)

		#
		# COMPLEX NUMBER PROJECTS
		#
		# If your project works with Complex numbers (STATE type = complex any type), please add it on the list below
		if(BUILD_COMPLEX_PROJECTS)

			#Example to work with scalar as complex
			#add_subdirectory(Helmholtz1D)
			add_subdirectory(TesteComplexo1D)
			add_subdirectory(StriplineFran)
			add_subdirectory(ValidacaoHCurlFran1)
			add_subdirectory(ValidacaoHCurlFran2)

			#add_subdirectory(HDiv_demo)

			#Examples used as tutorial for FEEC group (Prof Hugo)
			#add_subdirectory(PDElinear1D)
			#add_subdirectory(PDElinear2D)

		else(BUILD_COMPLEX_PROJECTS)
			#
			# REGULAR PROJECTS
			#
			# If your project DOES NOT use Complex numbers, please add it on the list below
			#
			########################## Nao aplicavel para Windows
			if (NOT WIN32)
				#Projeto Teo Refinamento Paralelo
				add_subdirectory(RefThread)
				# Projeto Caju: precisamos consertar para rodar no windows! problema com SWX components
				#add_subdirectory(HydraulicFracturePropagation)
			endif()
			##########################
			#Testes Douglas
			add_subdirectory(dactests)

			#Using discontinuous approximation functions - masters thesis of Erick Slis
			add_subdirectory(descontinuo)


			#Testing generating ellipse from a vector of points by least square method
			add_subdirectory(EllipseGenerator)

			#Adaptivity based on mesh cloning (Edimar Cesar thesis)
			#add_subdirectory(Poisson_ArcTan)

			#add_subdirectory(projteste)

			add_subdirectory(MultipleLoad)

			# add_subdirectory(adaptivity_refactor)
						
			add_subdirectory(Cubo_do_Nathan)

			#add_subdirectory(PyramidHdivTests)
			
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
			#Add_subdirectory(Conductivity)
			add_subdirectory(ConsLaw)
			#add_subdirectory(Descontinuo2)

			add_subdirectory(SteamInjection)
			
			#projeto agnaldo 09-09-2011
			#add_subdirectory(DoubleProjection)
                        add_subdirectory(artigoPira)
			add_subdirectory(TrabalhosLNCC15)
			#projeto omar duran 16-12-2012
			add_subdirectory(PoroElasticAnalysis)
			add_subdirectory(Poisson2DGID)
			add_subdirectory(OilWaterSystem)
			#add_subdirectory(DarcyflowAxisymmetricHdiv)
			#add_subdirectory(DarcyflowHdiv3D)
			add_subdirectory(MonophasicWell)
			if(USING_FAD)
				add_subdirectory(GradOfX)
			endif()
			#collaborative project (NH GH PD, SG, OD) 18-04-2015
			#add_subdirectory(iRMS)

			add_subdirectory(Elasticity2D)
			add_subdirectory(Gridgeneration)
			#add_subdirectory(Poisson2DLabyrinth)
			#add_subdirectory(Andrade)
			add_subdirectory(MatrixMarket)
			add_subdirectory(Frac1DHdiv)
			add_subdirectory(SubStruct)
			add_subdirectory(MultiphysicsSimulation)
			add_subdirectory(ProblemGirkmann)

			add_subdirectory(multplaca)

			#Projeto Lua 2012
			add_subdirectory(Lua01)

			#teste nathan
			add_subdirectory(PlasticityPV)

			#Projeto Denise
			#add_subdirectory(HdivPressureTeste)
			add_subdirectory(HpAdaptivity)
			add_subdirectory(Hdiv3DCurved)
			add_subdirectory(Hdiv3D)
			add_subdirectory(AcoplamentoH1Hdiv)
			add_subdirectory(Validacao)

			#add_subdirectory(prodindex)
			#Projeto Agnaldo
			add_subdirectory(Testes_Agnaldo)
			add_subdirectory(ReducedSpace)
			add_subdirectory(ReducedSpacePlasticity)
			add_subdirectory(HydraulicFracturePropagation)
			add_subdirectory(StabilizedHybridElasticity)
			add_subdirectory(GlobalPostPro)
			add_subdirectory(HDiv)
			add_subdirectory(HDiv_demo)

			# validacaoFranHCurl
			add_subdirectory(ValidacaoHCurlFran1)

			#Projetos de Thiago Dias dos Santos
			#add_subdirectory(SantosProjects)

			#GPU project
			add_subdirectory(GPU)

			#Projeto Pablo
			add_subdirectory(CreepStokes)
			
			add_subdirectory(Bima1)
#			add_subdirectory(EstabSUPG)
#			add_subdirectory(MetodosEstab)
			#
			# Colaboration with Laurent Stainier

			add_subdirectory(MatLib)

			add_subdirectory(MHM)

			if(USING_MKL)
				add_subdirectory(PardisoSample)
			endif(USING_MKL)

		endif(BUILD_COMPLEX_PROJECTS)
		
	endif(BUILD_PROJECTS)
