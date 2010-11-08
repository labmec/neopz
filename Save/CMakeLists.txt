
file(GLOB headers *.h)
install(FILES ${headers} DESTINATION include)

file(GLOB sources *.cpp)

########### next target ###############

set(save_STAT_SRCS
	${headers}
	${sources}
	PARENT_SCOPE
	)


########### next target ###############


#original Makefile.am contents follow:

##
## NeoPZ [ directory ] makefile
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
#savedir = $(prefix)/Save
#
#saveincdir = $(includedir)/Save
#saveinc_HEADERS= $(top_srcdir)/Save/*.h
#
#
#noinst_LIBRARIES = libSave.a
#
#libSave_a_SOURCES = pzsave.cpp 	pzfilebuffer.cpp pzbfilestream.cpp
#
### --| NeoPZ |------------------------------------------------------------------
