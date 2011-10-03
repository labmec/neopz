########### next target ###############

file(GLOB header *.h)

set(main_SRCS conductivity.cpp pzvoidflux.cpp)

set(flux_SRCS pzvoidflux.cpp voidflux_test.cpp TPZConductivityProblem.cpp)

set(conduct_SRCS pzvoidflux.cpp TPZConductivityProblem.cpp TPZConductivityMain.cpp)

add_executable(conductivity ${main_SRCS} ${header})

if(USING_BOOST)
add_executable(VoidFluxTest ${flux_SRCS} ${header})
endif()

add_executable(VoidFlux ${conduct_SRCS} ${header})

target_link_libraries(conductivity pz pthread)

if(USING_BOOST)
target_link_libraries(VoidFluxTest pz pthread boost_unit_test_framework)
endif()

target_link_libraries(VoidFlux pz pthread)

