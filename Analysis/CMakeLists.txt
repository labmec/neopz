
file(GLOB headers *.h)
install(FILES ${headers} DESTINATION include)

file(GLOB sources *.cpp)

########### next target ###############

set(pzanalysis_STAT_SRCS
	${headers}
	${sources}
	PARENT_SCOPE
	)




#original Makefile.am contents follow:

##
## NeoPZ Analysis makefile
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
#analysisdir = $(prefix)/Analysis
#
#analysisincdir = $(includedir)/Analysis
#analysisinc_HEADERS= $(top_srcdir)/Analysis/*.h
#
#noinst_LIBRARIES = libpzanalysis.a
#
#libpzanalysis_a_SOURCES = TPZNLMultGridAnalysis.cpp TPZNLMultGridAnalysis.h \
#	pzanalysis.cpp pzanalysiserror.cpp pzblackoilanalysis.cpp pzeuleranalysis.cpp \
#	pzmganalysis.cpp pznonlinanalysis.cpp pzsmanal.cpp pzsmfrontalanal.cpp \
#	pztransientanalysis.cpp
#
### --| NeoPZ |------------------------------------------------------------------
#
