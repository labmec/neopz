# @file neopz/Projects/CMakeList.txt  -- First CMake file for Projects directory 


if(BUILD_PROJECTS)
	
	#add_subdirectory(Cubo_do_Nathan)
	
	#add_subdirectory(HyperElasticityFAD)

	#Increment singular shape function: SingularShapeFunction
	#add_subdirectory(SingularShapeFunction)
	#add_subdirectory(angletest)
	#add_subdirectory(checkconv)

	#Orthogonalization process of Gram-Schmidt
	#add_subdirectory(matrix)

	#Inserting material and boundary conditions by function
	#add_subdirectory(termica)
	#add_subdirectory(Steam)
	#add_subdirectory(QuadraticElements)
	#add_subdirectory(ThermalStress)
	#add_subdirectory(HDiv)
	#add_subdirectory(Conductivity)
	#add_subdirectory(ConsLaw)
	#add_subdirectory(Descontinuo2)
	#add_subdirectory(HydraulicFracturePropagation)

	#add_subdirectory(SteamInjection)
	
	#projeto agnaldo 09-09-2011
	#add_subdirectory(DoubleProjection1D)
	#add_subdirectory(SubStruct)
	#add_subdirectory(Perf01_SubStruct)
	#add_subdirectory(Perf02_SubStruct-Nathan)
        add_subdirectory(Perf02_SubStruct-Nathan-Timed)
	#add_subdirectory(Perf02_SubStruct-Merge)
	#add_subdirectory(MultiphysicsSimulation)
	#add_subdirectory(prodindex)
	
	
endif(BUILD_PROJECTS)
