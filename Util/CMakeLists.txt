
file(GLOB headers *.h)
install(FILES ${headers} DESTINATION include)

file(GLOB sources *.cpp)

########### next target ###############

set(util_STAT_SRCS
	${headers}
	${sources}
	PARENT_SCOPE
	)




#original Makefile.am contents follow:

##
## NeoPZ Util makefile
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
#
#utildir = $(prefix)/Util
#
#utilincdir = $(includedir)/Util
#utilinc_HEADERS= $(top_srcdir)/Util/*.h
#
#noinst_LIBRARIES = libpzutil.a
#
#libpzutil_a_SOURCES = TPZTimer.cpp pzfunction.cpp pzline.cpp pzlog.cpp \
#	pzmap.cpp pznumeric.cpp pzplane.cpp pzpolynomial.cpp pzstring.cpp \
#	tpzautopointer.cpp tpzpermutation.cpp
#
### --| NeoPZ |------------------------------------------------------------------
#noinst_HEADERS = pzfunction.h
