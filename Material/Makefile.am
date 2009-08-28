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

libpzmaterial_a_SOURCES = 	TPZDiffusionConsLaw.cpp 	TPZPlacaOrthotropic.cpp \
		swelling.cpp TPZConsLawTest.cpp TPZConservationLaw.cpp TPZDiffusionConsLaw.cpp \
	TPZEuler.cpp TPZEulerConsLaw.cpp TPZFVHybrid.cpp TPZLinearConvection.cpp \
	TPZMulticamadaOrtho.cpp eulerdif.cpp pzartdiff.cpp pzausmflux.cpp pzbctension.cpp \
	pzbiharmonic.cpp pzblackoil2p3d.cpp pzbndcond.cpp pzburger.cpp pzconslaw.cpp \
	pzcoupledtransportdarcy.cpp pzcoupledtransportdarcyBC.cpp pzdiscgal.cpp pzelasAXImat.cpp \
	pzelasmat.cpp pzelast3d.cpp pzeuler.cpp pzeulerconslaw.cpp pzgradientflux.cpp \
	pzincnskeps.cpp pzl2projection.cpp pzmat1dlin.cpp pzmat2dlin.cpp pzmaterial.cpp \
	pzmaterialdata.cpp pzmaterialid.cpp pzmathyperelastic.cpp pzmatorthotropic.cpp \
	pzmatplaca2.cpp pzmattest.cpp pzmattest3d.cpp pzmultplaca.cpp pznlmat1d.cpp \
	pznlmat1drotatedengstrain.cpp pznonlinbiharmonic.cpp pznonlinearpoisson3d.cpp pzplaca.cpp \
	pzpoisson3d.cpp pzpoisson3dreferred.cpp pzspacetimerichardseq.cpp pzthermicelast3d.cpp \
	pztransientmat.cpp tpzmultcamada.cpp tpzoutofrange.cpp

## --| NeoPZ |------------------------------------------------------------------
noinst_HEADERS = pzausmflux.h pzburger.h pzeuler.h pzgradientflux.h \
	pzl2projection.h pzmaterialdata.h pzspacetimerichardseq.h

