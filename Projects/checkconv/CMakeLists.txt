
########### next target ###############

set(main_SRCS convhyper.cpp main.cpp)

add_executable(checkconv ${main_SRCS})

target_link_libraries(checkconv pz pthread)


########### install files ###############




#original Makefile.am contents follow:

##
## NeoPZ checkconv makefile
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
#main_SOURCES = convhyper.cpp main.cpp
#
### --| NeoPZ |------------------------------------------------------------------
