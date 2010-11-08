
include_directories(${KDE4_INCLUDES} ${KDE4_INCLUDE_DIR} ${QT_INCLUDES} )


########### next target ###############

set(RefPat_SRCS main.cpp)

kde4_add_executable(RefPat ${RefPat_SRCS})

target_link_libraries(RefPat ${KDE4_KDECORE_LIBS} pz metis)

install(TARGETS RefPat ${INSTALL_TARGETS_DEFAULT_ARGS})


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
#bin_PROGRAMS = RefPat
#
#RefPat_SOURCES = main.cpp
#
#RefPat_LDADD =  $(top_builddir)/lib/libpz.a @METISDIR@/lib/libmetis.a 
#
### --| NeoPZ |------------------------------------------------------------------
