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

libpzmesh_a_SOURCES = pzelctemp.cpp \
	pzcheckgeom.cpp \
	pzcheckmesh.cpp \
	pzcheckrestraint.cpp \
	pzcmesh.cpp \
        pzflowcmesh.cpp \
	pzcompel.cpp \
	pzconnect.cpp \
	pzelc1d.cpp \
	pzelcc3d.cpp \
	pzelcpi3d.cpp \
	pzelcpoint.cpp \
	pzelcpr3d.cpp \
	pzelcq2d.cpp \
	pzelct2d.cpp \
	pzelct3d.cpp \
	pzelg1d.cpp \
	pzelgc3d.cpp \
	pzelgpi3d.cpp \
	pzelgpoint.cpp \
	pzelgpr3d.cpp \
	pzelgq2d.cpp \
	pzelgt2d.cpp \
	pzelgt3d.cpp \
	pzelmat.cpp \
	pzgeoelbc.cpp \
	pzgeoel.cpp \
	pzgeoelside.cpp \
	pzgeoelrefless.cpp \
	pzgmesh.cpp \
	pzgnode.cpp \
	pzintel.cpp \
	pzstackmesh.cpp \
	pzsubcmesh.cpp \
	pzvecmesh.cpp \
	TPZCompElDisc.cpp \
	TPZInterfaceEl.cpp \
	TPZAgglomerateEl.cpp

## --| NeoPZ |------------------------------------------------------------------
