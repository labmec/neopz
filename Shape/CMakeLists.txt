
file(GLOB headers *.h)
install(FILES ${headers} DESTINATION include)

file(GLOB sources *.cpp)

########### next target ###############

set(pzshape_STAT_SRCS
	${headers}
	${sources}
	PARENT_SCOPE
	)



#original Makefile.am contents follow:

##
## NeoPZ Shape makefile
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
#shapedir = $(prefix)/Shape
#
#shapeincdir = $(includedir)/Shape
#shapeinc_HEADERS= $(top_srcdir)/Shape/*.h
#
#noinst_LIBRARIES = libpzshape.a
#
#libpzshape_a_SOURCES = pzshapetriang.cpp pzshapecube.cpp pzshapelinear.cpp pzshapepiram.cpp pzshapeprism.cpp pzshapequad.cpp pzshapetetra.cpp  TPZShapeDisc.cpp 
#
### --| NeoPZ |------------------------------------------------------------------
