
include_directories(${KDE4_INCLUDES} ${KDE4_INCLUDE_DIR} ${QT_INCLUDES} )

add_subdirectory(Class_5)
add_subdirectory(Class_1)
add_subdirectory(Class_2)
add_subdirectory(Class_3)
add_subdirectory(Class_4)


########### install files ###############




#original Makefile.am contents follow:

##
## NeoPZ Projects makefile
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
#SUBDIRS = Class_5 Class_1 Class_2 Class_3 Class_4
#
### --| NeoPZ |------------------------------------------------------------------
