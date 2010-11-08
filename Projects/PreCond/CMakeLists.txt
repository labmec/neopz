
include_directories(${KDE4_INCLUDES} ${KDE4_INCLUDE_DIR} ${QT_INCLUDES} )


########### next target ###############

set(PreCond_SRCS PreCond.cpp)

kde4_add_executable(PreCond ${PreCond_SRCS})

target_link_libraries(PreCond ${KDE4_KDECORE_LIBS})

install(TARGETS PreCond ${INSTALL_TARGETS_DEFAULT_ARGS})


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
#bin_PROGRAMS = PreCond
#
#PreCond_SOURCES = PreCond.cpp 
#
#
### --| NeoPZ |------------------------------------------------------------------
