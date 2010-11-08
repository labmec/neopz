
include_directories(${KDE4_INCLUDES} ${KDE4_INCLUDE_DIR} ${QT_INCLUDES} )


########### next target ###############

set(jornadas_SRCS jornadas.cpp)

kde4_add_executable(jornadas ${jornadas_SRCS})

target_link_libraries(jornadas ${KDE4_KDECORE_LIBS})

install(TARGETS jornadas ${INSTALL_TARGETS_DEFAULT_ARGS})


########### install files ###############




#original Makefile.am contents follow:

##
## NeoPZ multplaca makefile
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
#bin_PROGRAMS = jornadas
#
#jornadas_SOURCES = jornadas.cpp
#
### --| NeoPZ |------------------------------------------------------------------
