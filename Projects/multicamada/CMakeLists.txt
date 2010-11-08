
include_directories(${KDE4_INCLUDES} ${KDE4_INCLUDE_DIR} ${QT_INCLUDES} )


########### next target ###############

set(main_SRCS multicamada.cpp)

kde4_add_executable(main ${main_SRCS})

target_link_libraries(main ${KDE4_KDECORE_LIBS} pz metis)

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
#main_SOURCES = multicamada.cpp
#
#
#main_LDADD = $(top_builddir)/lib/libpz.a @METISDIR@/lib/libmetis.a
#
### --| NeoPZ |------------------------------------------------------------------
