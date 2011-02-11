
########### next target ###############

set(main_SRCS mainshape.cpp)

add_executable(TestShape ${main_SRCS})

target_link_libraries(TestShape pz pthread)

enable_testing()

add_test(TestShapeRuns TestShape)


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
#main_SOURCES = mainshape.cpp
#
### --| NeoPZ |------------------------------------------------------------------
#main_DEPENDENCIES = $(top_builddir)/lib/libpz.a
