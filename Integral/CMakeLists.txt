
file(GLOB headers *.h)
install(FILES ${headers} DESTINATION include)

file(GLOB sources *.cpp)

########### next target ###############

set(pzintegral_STAT_SRCS
	${headers}
	${sources}
	PARENT_SCOPE
	)



#original Makefile.am contents follow:

##
## NeoPZ Integral makefile
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
#integraldir = $(prefix)/Integral
#
#integralincdir = $(includedir)/Integral
#integralinc_HEADERS= $(top_srcdir)/Integral/*.h
#
#
#noinst_LIBRARIES = libpzintegral.a
#
#libpzintegral_a_SOURCES = pzquad.cpp tpzint1point.cpp tpzintpoints.cpp \
#    tpzgaussrule.cpp tpzgausslobattorule.cpp tpzintrule.cpp tpzintrulelist.cpp \
#    tpzintrulep3d.cpp tpzintrulet.cpp tpzintrulet3d.cpp tpzprinteg.cpp
#
### --| NeoPZ |------------------------------------------------------------------
#noinst_HEADERS = tpzint1point.h tpzintpoint.h tpzintpoints.h tpzintrule.h \
#	tpzintrulelist.h tpzintrulep3d.h tpzintrulet.h tpzintrulet3d.h tpzprinteg.h
