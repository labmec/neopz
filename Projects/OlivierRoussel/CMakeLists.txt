
include_directories(${KDE4_INCLUDES} ${KDE4_INCLUDE_DIR} ${QT_INCLUDES} )


########### next target ###############

set(lara_SRCS
    MultiResMesh.cpp
    adapt.cpp
    main.cpp
    malhas.cpp
    pzexplfinvolanal.cpp
    pzexplfinvolanal.h
    tools.cpp)

kde4_add_executable(lara ${lara_SRCS})

target_link_libraries(lara ${KDE4_KDECORE_LIBS})

install(TARGETS lara ${INSTALL_TARGETS_DEFAULT_ARGS})


########### install files ###############




#original Makefile.am contents follow:

#DISTCLEANFILES = @PZ_DISTCLEANFILES@
#CLEANFILES     = @PZ_CLEANFILES@
#
#-include $(top_builddir)/@PZ_MK@
#
#bin_PROGRAMS = lara 
#
#
#
### --| NeoPZ |------------------------------------------------------------------
#
#lara_SOURCES = MultiResMesh.cpp adapt.cpp main.cpp malhas.cpp \
#	pzexplfinvolanal.cpp pzexplfinvolanal.h tools.cpp
#
#lara_DEPENDENCIES = $(top_builddir)/lib/libpz.a
#
## noinst_HEADERS = TPZFakeFunction.h
## _SOURCES = adapt.h
#noinst_HEADERS = MultiResMesh.h
