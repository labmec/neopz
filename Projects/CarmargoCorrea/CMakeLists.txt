
include_directories(${KDE4_INCLUDES} ${KDE4_INCLUDE_DIR} ${QT_INCLUDES} )


########### next target ###############

set(cccc_SRCS ctriangletests.cpp)

kde4_add_executable(cccc ${cccc_SRCS})

target_link_libraries(cccc ${KDE4_KDECORE_LIBS})

install(TARGETS cccc ${INSTALL_TARGETS_DEFAULT_ARGS})


########### install files ###############




#original Makefile.am contents follow:

#
#DISTCLEANFILES = @PZ_DISTCLEANFILES@
#CLEANFILES     = @PZ_CLEANFILES@
#
#-include $(top_builddir)/@PZ_MK@ 
#
#INCLUDES = 
#METASOURCES = AUTO
#bin_PROGRAMS = cccc
#cccc_SOURCES = ctriangletests.cpp
#
#cccc_DEPENDENCIES = $(top_builddir)/lib/libpz.a
