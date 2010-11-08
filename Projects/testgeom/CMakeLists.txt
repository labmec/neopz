########### next target ###############

set(main_SRCS pzgbmesh.cpp geom.cpp)

add_executable(TestGeom ${main_SRCS})

configure_file(quad_st_800_1R_16X.gri quad_st_800_1R_16X.gri COPYONLY)

target_link_libraries(TestGeom pz)

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
#main_SOURCES = pzgbmesh.cpp geom.cpp 
#
#
### --| NeoPZ |------------------------------------------------------------------
#main_DEPENDENCIES = $(top_builddir)/lib/libpz.a
