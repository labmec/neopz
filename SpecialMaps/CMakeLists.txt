
file(GLOB headers *.h)
install(FILES ${headers} DESTINATION include)

file(GLOB sources *.cpp)

########### next target ###############

set(pzspecialmaps_STAT_SRCS
	${headers}
	${sources}
	PARENT_SCOPE
	)





#original Makefile.am contents follow:

##
## NeoPZ SpecialMaps makefile
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
#specialmapsdir = $(prefix)/SpecialMaps
#
#specialmapsincdir = $(includedir)/SpecialMaps
#specialmapsinc_HEADERS= $(top_srcdir)/SpecialMaps/*.h
#
#
#noinst_LIBRARIES = libpzspecialmaps.a
#
#libpzspecialmaps_a_SOURCES = convtest.cpp tpzarc3d.cpp tpzblendnaca.cpp \
#	tpzchangeel.cpp tpzcurvedtriangle.cpp tpzellipse.cpp tpzmathtools.cpp \
#	tpzquadraticline.cpp tpzquadraticquad.cpp tpzquadratictetra.cpp tpzquadratictrig.cpp
#
### --| NeoPZ |------------------------------------------------------------------
#noinst_HEADERS = convtest.h tpzarc3d.h tpzblendnaca.h tpzchangeel.h \
#	tpzcurvedtriangle.h tpzellipse.h tpzmathtools.h tpzquadraticline.h tpzquadraticquad.h \
#	tpzquadratictetra.h tpzquadratictrig.h
