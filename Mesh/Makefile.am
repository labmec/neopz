#
# NeoPZ Mesh makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id#
#

DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_builddir)/@PZ_MK@

meshdir = $(prefix)/Mesh
meshincdir = $(includedir)/Mesh
meshinc_HEADERS= $(top_srcdir)/Mesh/*.h

noinst_LIBRARIES = libpzmesh.a

libpzmesh_a_SOURCES = pzinterpolationspace.cpp pzcompel.cpp pzelctemp.cpp 	pzmeshid.cpp	 \
		TPZCompElDisc.cpp 	TPZInterfaceEl.cpp 	TPZAgglomerateEl.cpp 	pzgmesh.cpp 	pzcheckgeom.cpp \
			pzcheckmesh.cpp 	pzcheckrestraint.cpp 	pzcmesh.cpp	 tpzcompmeshreferred.cpp	 \
		pzflowcmesh.cpp 	pzconnect.cpp	pzelmat.cpp 	pzgeoelbc.cpp 	pzgeoel.cpp 	pzgeoelside.cpp \
			pzgeoelrefless.cpp 	pzgnode.cpp 	pzintel.cpp 	pzstackmesh.cpp 	pzvecmesh.cpp \
			tpzgeoelrefpattern.cpp tpzagglomeratemesh.cpp pzreferredcompel.cpp tpzgeoelmapped.cpp pzsubcmesh.cpp

## --| NeoPZ |------------------------------------------------------------------
noinst_HEADERS = pzgeoelrefless.h.h pzinterpolationspace.h tpzgeoelmapped.h \
	tpzgeoelrefpattern.h.h
