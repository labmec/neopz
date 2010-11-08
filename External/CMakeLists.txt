
add_subdirectory(sloan)


file(GLOB headers *.h)
install(FILES ${headers} DESTINATION include)

file(GLOB sources *.cpp)

########### next target ###############

set(pzexternal_STAT_SRCS
	${pzexternal_STAT_SRCS}
	${headers}
	${sources}
	PARENT_SCOPE
	)



#original Makefile.am contents follow:

##
## NeoPZ External makefile
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
#externaldir = $(prefix)/External
#
#externalincdir = $(includedir)/External
#externalinc_HEADERS= $(top_srcdir)/External/*.h
#
#
#noinst_LIBRARIES = libpzexternal.a
#
#libpzexternal_a_SOURCES = TPZBoostGraph.cpp TPZBoostGraph.h pzmetis.cpp \
#	pzrenumbering.cpp pzsloan.cpp tpznodesetcompute.cpp
#
### --| NeoPZ |------------------------------------------------------------------
#SUBDIRS = sloan
