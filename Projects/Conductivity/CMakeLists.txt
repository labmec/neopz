########### next target ###############

file(GLOB header *.h)

set(main_SRCS conductivity.cpp pzvoidflux.cpp)

set(flux_SRCS pzvoidflux.cpp voidflux_test.cpp TPZConductivityProblem.cpp)

set(conduct_SRCS pzvoidflux.cpp TPZConductivityProblem.cpp TPZConductivityMain.cpp)

# First executable: main function in conductivity.cpp
add_executable(conductivity ${main_SRCS} ${header})
target_link_libraries(conductivity pz pthread)

# Second executable: main function in TPZConductivityMain.cpp
add_executable(conductivity_VoidFlux ${conduct_SRCS} ${header})
target_link_libraries(conductivity_VoidFlux pz pthread)

# Third executable: main function in voidflux_test.cpp - It is a Unit Test project????
add_executable(conductivity_VoidFluxTest ${flux_SRCS} ${header})
if(USING_BOOST)
	target_link_libraries(conductivity_VoidFluxTest pz pthread boost_unit_test_framework)
else()
	target_link_libraries(conductivity_VoidFluxTest pz pthread)
endif()


