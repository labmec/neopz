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

libpzmaterial_a_SOURCES = eulerdif.cpp \
	pzbndcond.cpp \
	pzelasmat.cpp \
	pzmat1dlin.cpp \
	pzmat2dlin.cpp \
	pzmaterial.cpp \
	pzmathyperelastic.cpp \
	pzmatorthotropic.cpp \
	pzmatplaca2.cpp \
	pzmatplaca.cpp \
	pzmattest3d.cpp \
	pzmattest.cpp \
	pzmultplaca.cpp \
	pzplaca.cpp \
	pzpoisson3d.cpp \
	TPZEuler.cpp \
	TPZEulerConsLaw.cpp \
	TPZDiffusionConsLaw.cpp \
	TPZLinearConvection.cpp \
	tpzmultcamada.cpp \
	TPZMulticamadaOrtho.cpp \
	TPZPlacaOrthotropic.cpp \
	TPZConservationLaw.cpp \
	TPZConsLawTest.cpp \
	TPZFVHybrid.cpp \
	pzartdiff.cpp \
	swelling.cpp \
	pzconslaw.cpp \
	pzeulerconslaw.cpp \
	pzartdiff.cpp \
	pzbctension.cpp \
	pzbiharmonic.cpp

## --| NeoPZ |------------------------------------------------------------------
