
file(GLOB headers *.h)
install(FILES ${headers} DESTINATION include)

file(GLOB sources *.cpp)

########### next target ###############

set(pzgeom_STAT_SRCS
	${headers}
	${sources}
	PARENT_SCOPE
	)


#original Makefile.am contents follow:

##
## NeoPZ Geom makefile
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
#geomdir = $(prefix)/Geom
#
#geomincdir = $(includedir)/Geom
#geominc_HEADERS= $(top_srcdir)/Geom/*.h
#
#
#noinst_LIBRARIES = libpzgeom.a
#
#libpzgeom_a_SOURCES = 	TPZGeoCube.cpp PrismExtend.cpp TPZGeoElement.cpp \
#	TPZGeoLinear.cpp pzgeopoint.cpp pzgeoprism.cpp pzgeopyramid.cpp pzgeoquad.cpp \
#	pzgeotetrahedra.cpp pzgeotriangle.cpp pznoderep.cpp
#
### --| NeoPZ |------------------------------------------------------------------
#
#noinst_HEADERS = TPZGeoElement.h.h pznoderep.h
