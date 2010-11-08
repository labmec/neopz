
include_directories(${KDE4_INCLUDES} ${KDE4_INCLUDE_DIR} ${QT_INCLUDES} )


########### next target ###############

set(]_STAT_SRCS ${CMAKE_CURRENT_BINARY_DIR}/dummy.cpp)

kde4_add_library(] STATIC ${]_STAT_SRCS})


########### install files ###############




#original Makefile.am contents follow:

##
## NeoPZ [ directory ] makefile
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
#noinst_LIBRARIES = lib[ name of the library ].a
#
#lib[ name of the library ]_a_SOURCES = 
#
### --| NeoPZ |------------------------------------------------------------------
