########### next target ###############

file(GLOB locheaders *.h)

file(GLOB locsources *.cpp)

set(main_SRCS ${locheaders} ${locsources})

add_executable(Steam ${main_SRCS})

#configure_file(quad_st_800_1R_16X.gri quad_st_800_1R_16X.gri COPYONLY)

target_link_libraries(Steam pz pthread)

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
