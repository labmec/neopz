#
# NeoPZ Material makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id#
#

DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_srcdir)/@PZ_MK@

materialdir = $(prefix)/Material

material_HEADERS = $(wildcard *.h)

noinst_LIBRARIES = libpzmaterial.a

libpzmaterial_a_SOURCES = eulerdif.cc \
	pzbndcond.cc \
	pzelasmat.cc \
	pzmat1dlin.cc \
	pzmat2dlin.cc \
	pzmaterial.cc \
	pzmathyperelastic.cc \
	pzmatorthotropic.cc \
	pzmatplaca2.cc \
	pzmatplaca.cc \
	pzmattest3d.cc \
	pzmattest.cc \
	pzmultplaca.cc \
	pzplaca.cc \
	pzpoisson3d.cc \
	TPZEuler.cc \
	TPZLinearConvection.cc \
	tpzmultcamada.cc \
	TPZMulticamadaOrtho.cc \
	TPZPlacaOrthotropic.cc 

## --| NeoPZ |------------------------------------------------------------------
