
include_directories(${KDE4_INCLUDES} ${KDE4_INCLUDE_DIR} ${QT_INCLUDES} )


########### next target ###############

set(main_SRCS
    main.cpp
    pzerror_ind.cpp
    tmbheidi.cpp
    tmbadapinterface.cpp
    tmbxubinshih.cpp)

kde4_add_executable(main ${main_SRCS})

target_link_libraries(main ${KDE4_KDECORE_LIBS})

install(TARGETS main ${INSTALL_TARGETS_DEFAULT_ARGS})


########### install files ###############




#original Makefile.am contents follow:

##
## NeoPZ error  makefile
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
#main_SOURCES = main.cpp pzerror_ind.cpp tmbheidi.cpp tmbadapinterface.cpp tmbxubinshih.cpp
#
### --| NeoPZ |------------------------------------------------------------------
