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

libpzmesh_a_SOURCES = pzelctemp.cc \
	pzcheckgeom.cc \
	pzcheckmesh.cc \
	pzcheckrestraint.cc \
	pzcmesh.cc \
        pzflowcmesh.cc \
	pzcompel.cc \
	pzconnect.cc \
	pzelc1d.cc \
	pzelcc3d.cc \
	pzelcpi3d.cc \
	pzelcpoint.cc \
	pzelcpr3d.cc \
	pzelcq2d.cc \
	pzelct2d.cc \
	pzelct3d.cc \
	pzelg1d.cc \
	pzelgc3d.cc \
	pzelgpi3d.cc \
	pzelgpoint.cc \
	pzelgpr3d.cc \
	pzelgq2d.cc \
	pzelgt2d.cc \
	pzelgt3d.cc \
	pzelmat.cc \
	pzgeoelbc.cc \
	pzgeoel.cc \
	pzgeoelside.cc \
	pzgmesh.cc \
	pzgnode.cc \
	pzintel.cc \
	pzstackmesh.cc \
	pzsubcmesh.cc \
	pzvecmesh.cc \
	TPZCompElDisc.cc \
	TPZInterfaceEl.cc \
	TPZAgglomerateEl.cc

## --| NeoPZ |------------------------------------------------------------------
