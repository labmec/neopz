
include_directories(${KDE4_INCLUDES} ${KDE4_INCLUDE_DIR} ${QT_INCLUDES} )


########### next target ###############

set(main_SRCS main.cc)

kde4_add_executable(main ${main_SRCS})

target_link_libraries(main ${KDE4_KDECORE_LIBS})

install(TARGETS main ${INSTALL_TARGETS_DEFAULT_ARGS})


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
#main_SOURCES = main.cc 
#
## /home/ic/luis/Documents/dxflib-2.0.1.3-1.src/lib/libdxf.a
#
### --| NeoPZ |------------------------------------------------------------------
