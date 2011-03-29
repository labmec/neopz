########### next target ###############

set(main_SRCS termi.cpp fast.cpp)

add_executable(termica ${main_SRCS})

target_link_libraries(termica pz pthread)


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
#main_SOURCES = termi.cpp fast.cpp
#
