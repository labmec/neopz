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

libpzmaterial_a_SOURCES = eulerdif.cpp 	pzbndcond.cpp 	pzelasmat.cpp 	pzmat1dlin.cpp 	pzmat2dlin.cpp 	pzmaterial.cpp 	pzmathyperelastic.cpp 	pzmatorthotropic.cpp 	pzmatplaca2.cpp 	pzmatplaca.cpp 	pzmattest3d.cpp 	pzmattest.cpp 	pzmultplaca.cpp 	pzplaca.cpp 	pzpoisson3d.cpp 	TPZEuler.cpp 	TPZEulerConsLaw.cpp 	TPZDiffusionConsLaw.cpp	TPZLinearConvection.cpp 	tpzmultcamada.cpp 	TPZMulticamadaOrtho.cpp 	TPZPlacaOrthotropic.cpp 	TPZConservationLaw.cpp 	TPZConsLawTest.cpp 	TPZFVHybrid.cpp 	pzartdiff.cpp 	swelling.cpp 	pzconslaw.cpp 	pzeulerconslaw.cpp 	pzartdiff.cpp 	pzbctension.cpp 	pzbiharmonic.cpp pzmaterialid.cpp tpzoutofrange.cpp

## --| NeoPZ |------------------------------------------------------------------
noinst_HEADERS = pzmaterialid.h tpzoutofrange.h
INCLUDES = -I$(top_srcdir)/Analysis -I$(top_srcdir)/Common -I$(top_srcdir)/Integral -I$(top_srcdir)/Matrix -I$(top_srcdir)/Mesh -I$(top_srcdir)/Multigrid -I$(top_srcdir)/Post -I$(top_srcdir)/Save -I$(top_srcdir)/Shape -I$(top_srcdir)/StrMatrix -I$(top_srcdir)/Util @FADINCLUDE_TFADET@ @FADINCLUDE_TFAD@ @FADINCLUDE_FAD@ @FADINCLUDE@
AM_CXXFLAGS = @PZ_CXXFLAGS@ @USING_METIS@ @USING_SLOAN@ @USING_FAD@ @USING_ATLAS@ @USING_BLAS@ @ATLASINCLUDE@ 
AM_CXXFLAGS = @PZ_CXXFLAGS@ @USING_METIS@ @USING_SLOAN@ @USING_FAD@ @USING_ATLAS@ @USING_BLAS@
