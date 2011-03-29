
########### next target ###############

set(main_SRCS main.cpp TExtFunction.cpp TDiscoFunction.cpp)

add_executable(discShape ${main_SRCS})

target_link_libraries(discShape pz pthread)

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
#main_SOURCES = main.cpp TExtFunction.cpp TDiscoFunction.cpp
#
#main_DEPENDENCIES = $(top_builddir)/lib/libpz.a
