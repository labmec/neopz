
########### next target ###############

set(main_SRCS viga3d.cpp)

kde4_add_executable(main ${main_SRCS})

target_link_libraries(main PZ)

install(TARGETS main ${INSTALL_TARGETS_DEFAULT_ARGS})


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
#main_SOURCES = viga3d.cpp
#
### --| NeoPZ |------------------------------------------------------------------
