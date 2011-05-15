########### next target ###############

set(main_SRCS conductivity.cpp pzvoidflux.cpp)

set(flux_SRCS pzvoidflux.cpp voidflux_test.cpp TPZConductivityProblem.cpp)

set(conduct_SRCS pzvoidflux.cpp TPZConductivityProblem.cpp TPZConductivityMain.cpp)

add_executable(conductivity ${main_SRCS})

if(USING_BOOST)
add_executable(VoidFluxTest ${flux_SRCS})
endif()

add_executable(VoidFlux ${conduct_SRCS})

target_link_libraries(conductivity pz pthread)

if(USING_BOOST)
target_link_libraries(VoidFluxTest pz pthread boost_unit_test_framework)
endif()

target_link_libraries(VoidFlux pz pthread)

########### install files ###############




#original Makefile.am contents follow:

##
## NeoPZ adaptivity makefile
##
## It will be processed by GNU automake to produce Makefile.in
##
## $Id#
##
#
#DISTCLEANFILES = @PZ_DISTCLEANFILES@
#CLEANFILES     = @PZ_CLEANFILES@
#
#-include $(top_builddir)/@PZ_MK@
#
#bin_PROGRAMS = main 
#
#
#
### --| NeoPZ |------------------------------------------------------------------
#
#main_SOURCES = termi.cpp fast.cpp
#
