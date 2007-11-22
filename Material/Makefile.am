#
# NeoPZ Material makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id#
#

DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_builddir)/@PZ_MK@

materialdir = $(prefix)/Material

materialincdir = $(includedir)/Material
materialinc_HEADERS= $(top_srcdir)/Material/*.h


noinst_LIBRARIES = libpzmaterial.a

libpzmaterial_a_SOURCES = pzmaterialdata.cpp pzdiscgal.cpp eulerdif.cpp \
	pzbndcond.cpp pzelasmat.cpp 	pzmat1dlin.cpp pzmat2dlin.cpp pzmaterial.cpp \
	pzmathyperelastic.cpp 	pzmatorthotropic.cpp pzmatplaca2.cpp pzmatplaca.cpp pzmattest3d.cpp \
	pzmattest.cpp 	pzmultplaca.cpp pzplaca.cpp pzpoisson3d.cpp TPZEuler.cpp \
	TPZEulerConsLaw.cpp 	TPZDiffusionConsLaw.cpp TPZLinearConvection.cpp tpzmultcamada.cpp \
	TPZMulticamadaOrtho.cpp 	TPZPlacaOrthotropic.cpp TPZConservationLaw.cpp TPZConsLawTest.cpp \
	TPZFVHybrid.cpp pzartdiff.cpp 	swelling.cpp pzconslaw.cpp pzeulerconslaw.cpp \
	pzartdiff.cpp pzbctension.cpp 	pzbiharmonic.cpp pznonlinbiharmonic.cpp \
	pzmaterialid.cpp tpzoutofrange.cpp 	pzincnskeps.cpp pzelast3d.cpp \
	pzcoupledtransportdarcy.cpp pzcoupledtransportdarcyBC.cpp 	pzpoisson3dreferred.cpp \
	pzthermicelast3d.cpp pznlmat1drotatedengstrain.cpp pznlmat1d.cpp 	pztransientmat.cpp \
	pznonlinearpoisson3d.cpp pzburger.cpp pzl2projection.cpp

## --| NeoPZ |------------------------------------------------------------------
noinst_HEADERS = pzmaterialdata.h pzl2projection.h pzburger.h
