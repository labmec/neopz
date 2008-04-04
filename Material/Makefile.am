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
		pzbiharmonic.cpp 	pzincnskeps.cpp 	pzmat1dlin.cpp 	pzmatorthotropic.cpp 	pzmultplaca.cpp \
		pzpoisson3dreferred.cpp 	pztransientmat.cpp 	swelling.cpp TPZConsLawTest.cpp \
	TPZConservationLaw.cpp TPZEuler.cpp TPZEulerConsLaw.cpp TPZFVHybrid.cpp \
	TPZLinearConvection.cpp TPZMulticamadaOrtho.cpp eulerdif.cpp pzartdiff.cpp pzbctension.cpp \
	pzbndcond.cpp pzburger.cpp pzconslaw.cpp pzcoupledtransportdarcy.cpp \
	pzcoupledtransportdarcyBC.cpp pzdiscgal.cpp pzelasmat.cpp pzelast3d.cpp pzeulerconslaw.cpp \
	pzl2projection.cpp pzmat2dlin.cpp pzmaterial.cpp pzmaterialdata.cpp pzmaterialid.cpp \
	pzmathyperelastic.cpp pzmatplaca2.cpp pzmattest.cpp pzmattest3d.cpp pznlmat1d.cpp \
	pznlmat1drotatedengstrain.cpp pznonlinbiharmonic.cpp pznonlinearpoisson3d.cpp pzplaca.cpp \
	pzpoisson3d.cpp pzspacetimerichardseq.cpp pzthermicelast3d.cpp tpzmultcamada.cpp \
	tpzoutofrange.cpp

## --| NeoPZ |------------------------------------------------------------------
noinst_HEADERS = pzburger.h pzl2projection.h pzmaterialdata.h \
	pzspacetimerichardseq.h
