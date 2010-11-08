
include_directories(${KDE4_INCLUDES} ${KDE4_INCLUDE_DIR} ${QT_INCLUDES} )


########### next target ###############

set(testeangular_SRCS testeangular.cpp)

kde4_add_executable(testeangular ${testeangular_SRCS})

target_link_libraries(testeangular ${KDE4_KDECORE_LIBS})

install(TARGETS testeangular ${INSTALL_TARGETS_DEFAULT_ARGS})


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
