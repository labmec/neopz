#
# NeoPZ Mesh makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id#
#

DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_srcdir)/@PZ_MK@

meshdir = $(prefix)/Mesh

mesh_HEADERS = $(wildcard *.h)

noinst_LIBRARIES = libpzmesh.a

libpzmesh_a_SOURCES = pzsubcmesh.cpp 	pzelctemp.cpp 	pzmeshid.cpp	 TPZCompElDisc.cpp 	TPZInterfaceEl.cpp 	TPZAgglomerateEl.cpp 	pzgmesh.cpp 	pzcheckgeom.cpp 	pzcheckmesh.cpp 	pzcheckrestraint.cpp 	pzcmesh.cpp         pzflowcmesh.cpp 	pzcompel.cpp 	pzconnect.cpp 	pzelc1d.cpp 	pzelcc3d.cpp 	pzelcpi3d.cpp 	pzelcpoint.cpp 	pzelcpr3d.cpp 	pzelcq2d.cpp 	pzelct2d.cpp 	pzelct3d.cpp 	pzelg1d.cpp 	pzelgc3d.cpp 	pzelgpi3d.cpp 	pzelgpoint.cpp 	pzelgpr3d.cpp 	pzelgq2d.cpp 	pzelgt2d.cpp 	pzelgt3d.cpp 	pzelmat.cpp 	pzgeoelbc.cpp 	pzgeoel.cpp 	pzgeoelside.cpp 	pzgeoelrefless.cpp 	pzgnode.cpp 	pzintel.cpp 	pzstackmesh.cpp 	pzvecmesh.cpp 	tpzgeoelrefpattern.cpp  tpzagglomeratemesh.cpp

## --| NeoPZ |------------------------------------------------------------------
INCLUDES = -I$(top_srcdir)/Analysis -I$(top_srcdir)/Common -I$(top_srcdir)/Cosys -I$(top_srcdir)/External -I$(top_srcdir)/Frontal -I$(top_srcdir)/Geom -I$(top_srcdir)/Integral -I$(top_srcdir)/LinearSolvers -I$(top_srcdir)/Material -I$(top_srcdir)/Matrix -I$(top_srcdir)/Mesh -I$(top_srcdir)/Multigrid -I$(top_srcdir)/Post -I$(top_srcdir)/Refine -I$(top_srcdir)/Save -I$(top_srcdir)/Shape -I$(top_srcdir)/StrMatrix -I$(top_srcdir)/Util @FADINCLUDE_TFADET@ @FADINCLUDE_TFAD@ @FADINCLUDE_FAD@ @FADINCLUDE@
AM_CXXFLAGS = @PZ_CXXFLAGS@ @USING_METIS@ @USING_SLOAN@ @USING_FAD@ @USING_ATLAS@ @USING_BLAS@
