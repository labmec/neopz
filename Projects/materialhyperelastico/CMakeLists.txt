
include_directories(${KDE4_INCLUDES} ${KDE4_INCLUDE_DIR} ${QT_INCLUDES} )


########### next target ###############

set(hyperelastic_SRCS hyperelastic.cpp)

kde4_add_executable(hyperelastic ${hyperelastic_SRCS})

target_link_libraries(hyperelastic ${KDE4_KDECORE_LIBS})

install(TARGETS hyperelastic ${INSTALL_TARGETS_DEFAULT_ARGS})


########### install files ###############




#original Makefile.am contents follow:

##
## NeoPZ materialhyperelastico makefile
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
#bin_PROGRAMS = hyperelastic
#
#hyperelastic_SOURCES = hyperelastic.cpp
#
#
### --| NeoPZ |------------------------------------------------------------------
