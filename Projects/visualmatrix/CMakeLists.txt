
########### next target ###############

set(main_SRCS main.cpp)

add_executable(visualmatrix ${main_SRCS})

target_link_libraries(visualmatrix pz pthread)


########### install files ###############




#original Makefile.am contents follow:

##
## NeoPZ visualmatrix makefile
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
#main_SOURCES = main.cpp
#
### --| NeoPZ |------------------------------------------------------------------
