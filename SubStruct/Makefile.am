#
# NeoPZ SubStruct makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id#
#

DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_builddir)/@PZ_MK@

matrixdir = $(prefix)/SubStruct

matrixincdir = $(includedir)/SubStruct
matrixinc_HEADERS= $(top_srcdir)/SubStruct/*.h

noinst_LIBRARIES = libpzsubstruct.a

libpzsubstruct_a_SOURCES = pzdohrstructmatrix.cpp tpzdohrsubstruct.cpp tpzdohrsubstructCondense.cpp tpzdohrassemblelist.cpp tpzgensubstruct.cpp tpzdohrassembly.cpp tpzmatredstructmatrix.cpp tpzdohrmatrix.cpp tpzpairstructmatrix.cpp tpzdohrprecond.cpp TPZTimeTemp.cpp TPZfTime.cpp PrintTime.cpp TPZDataBase.cpp


## --| NeoPZ |------------------------------------------------------------------
