########### next target ###############
file(GLOB sources *.cpp)

set(main_SRCS ${sources})

add_executable(PerkinsGonzalez ${main_SRCS})

target_link_libraries(PerkinsGonzalez pz pthread)

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
