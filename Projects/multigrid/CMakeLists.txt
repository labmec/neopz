
include_directories(${KDE4_INCLUDES} ${KDE4_INCLUDE_DIR} ${QT_INCLUDES} )


########### next target ###############

set(multigrid_SRCS main.cpp)

kde4_add_executable(multigrid ${multigrid_SRCS})

target_link_libraries(multigrid ${KDE4_KDECORE_LIBS} pthread)

install(TARGETS multigrid ${INSTALL_TARGETS_DEFAULT_ARGS})


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
#bin_PROGRAMS = multigrid
#
#multigrid_SOURCES = main.cpp
#
#multigrid_LDADD = $(PZ_LIBDIR) $(PZ_LIB) -lpthread 
#
### --| NeoPZ |------------------------------------------------------------------
