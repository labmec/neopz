
set(testeangular_SRCS testeangular.cpp)

add_executable(testeangular ${testeangular_SRCS})

target_link_libraries(testeangular pz pthread)


########### install files ###############




#original Makefile.am contents follow:

##
## NeoPZ angletest makefile
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
#bin_PROGRAMS = testeangular
#
#testeangular_SOURCES = testeangular.cpp
#
### --| NeoPZ |------------------------------------------------------------------
