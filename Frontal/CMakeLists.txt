file(GLOB headers *.h)
install(FILES ${headers} DESTINATION include)

file(GLOB sources *.cpp)

########### next target ###############

set(pzfrontal_STAT_SRCS
	${headers}
	${sources}
	PARENT_SCOPE
	)




#original Makefile.am contents follow:

##
## NeoPZ Frontal makefile
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
#frontaldir = $(prefix)/Frontal
#
#frontalincdir = $(includedir)/Frontal
#frontalinc_HEADERS= $(top_srcdir)/Frontal/*.h
#
#
#noinst_LIBRARIES = libpzfrontal.a
#
#libpzfrontal_a_SOURCES = pzvecfrontal.cpp tpzeqnarray.cpp TPZFileEqnStorage.cpp TPZFront.cpp TPZFrontMatrix.cpp TPZFrontNonSym.cpp TPZFrontSym.cpp TPZParFrontMatrix.cpp TPZStackEqnStorage.cpp 
#
### --| NeoPZ |------------------------------------------------------------------
