########### next target ###############

file(GLOB header *.h)

set(main_SRCS conductivity.cpp pzvoidflux.cpp)

set(flux_SRCS pzvoidflux.cpp voidflux_test.cpp TPZConductivityProblem.cpp)

set(conduct_SRCS pzvoidflux.cpp TPZConductivityProblem.cpp TPZConductivityMain.cpp)

# First executable: main function in conductivity.cpp
add_executable(conductivity ${main_SRCS} ${header})
target_link_libraries(conductivity pz)

# Second executable: main function in TPZConductivityMain.cpp
add_executable(conductivity_VoidFlux ${conduct_SRCS} ${header})
target_link_libraries(conductivity_VoidFlux pz)

# Third executable: main function in voidflux_test.cpp - It is a Unit Test project????
if(BUILD_UNITTESTING)
	add_executable(conductivity_VoidFluxTest ${flux_SRCS} ${header})
	if(USING_BOOST)
		IF (WIN32)
			target_link_libraries(conductivity_VoidFluxTest pz ${Boost_LIBRARIES})
		ELSE()
			target_link_libraries(conductivity_VoidFluxTest pz)
		ENDIF()
	else()
		target_link_libraries(conductivity_VoidFluxTest pz)
	endif()
endif()
include (CTestTestFile.cmake)


