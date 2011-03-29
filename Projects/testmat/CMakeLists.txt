
########### next target ###############

set(main_SRCS mainmat.cpp)

add_executable(testmat ${main_SRCS})

target_link_libraries(testmat pz pthread)


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
#main_SOURCES = mainmat.cpp
#
#
### --| NeoPZ |------------------------------------------------------------------
