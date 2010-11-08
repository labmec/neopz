

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION include)

file(GLOB sources *.cpp)

########### next target ###############

set(common_STAT_SRCS
	${headers}
	${sources}
	PARENT_SCOPE
	)




#original Makefile.am contents follow:

##
## NeoPZ Common makefile
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
#commondir = $(prefix)/Common
#
#commonincdir = $(includedir)/Common
#commoninc_HEADERS= $(top_srcdir)/Common/*.h
#
#noinst_LIBRARIES = libCommon.a
#
#libCommon_a_SOURCES = pzreal.cpp
#
#
### --| NeoPZ |------------------------------------------------------------------
