#
# NeoPZ lib makefile
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id#
#

DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

-include $(top_builddir)/@PZ_MK@

lib_LIBRARIES = libpz.a

libpz_a_SOURCES =

libpz_a_LIBADD = ../External/pzsloan.o ../Util/TPZTimer.o \
	../Matrix/tpzsparseblockdiagonal.o ../External/tpznodesetcompute.o ../Common/pzreal.o ../Util/pzmap.o \
	../Util/pzadmchunkgen.o ../Util/pzchunkgen.o ../Util/pzavlmapgen.o ../Util/pzstackgen.o \
	../Util/pzmanvectorgen.o ../Util/pzstring.o ../Util/pzvecgen.o ../Refine/pzrefpoint.o \
	../Refine/pzrefprism.o ../Refine/pzrefpyram.o ../Refine/pzrefquad.o \
	../Refine/pzreftetrahedra.o ../Refine/pzreftriangle.o ../Refine/TPZRefCube.o \
	../Refine/TPZRefLinear.o ../Refine/TPZRefPattern.o ../Mesh/pzcmesh.o ../Mesh/pzflowcmesh.o \
	../Mesh/pzcompel.o ../Mesh/pzconnect.o ../Mesh/pzelc1d.o ../Mesh/pzcheckgeom.o \
	../Mesh/pzcheckmesh.o ../Mesh/pzcheckrestraint.o ../Mesh/pzelcc3d.o ../Mesh/pzelcpi3d.o \
	../Mesh/pzelcpoint.o ../Mesh/pzelcpr3d.o ../Mesh/pzelcq2d.o ../Mesh/pzelct2d.o \
	../Mesh/pzelct3d.o ../Mesh/pzelg1d.o ../Mesh/pzelgc3d.o ../Mesh/pzelgpi3d.o \
	../Mesh/pzelgpoint.o ../Mesh/pzelgpr3d.o ../Mesh/pzelgq2d.o ../Mesh/pzelgt2d.o \
	../Mesh/pzelgt3d.o ../Mesh/pzelmat.o ../Mesh/pzgeoelbc.o ../Mesh/pzgeoel.o \
	../Mesh/pzgeoelside.o ../Mesh/pzgeoelrefless.o ../Mesh/pzgmesh.o ../Mesh/pzgnode.o \
	../Mesh/pzintel.o ../Mesh/pzstackmesh.o ../Mesh/pzsubcmesh.o ../Mesh/pzvecmesh.o \
	../Mesh/TPZCompElDisc.o ../Mesh/TPZInterfaceEl.o ../Mesh/TPZAgglomerateEl.o ../Mesh/pzelctemp.o \
	../Mesh/tpzgeoelrefpattern.o ../Mesh/pzmeshid.o ../Shape/pzshapecube.o ../Shape/pzshapelinear.o \
	../Shape/pzshapepiram.o ../Shape/pzshapeprism.o ../Shape/pzshapequad.o ../Shape/pzshapetriang.o \
	../Shape/pzshapetetra.o ../Shape/TPZShapeDisc.o ../Matrix/pzlink.o ../Matrix/pzblock.o \
	../Matrix/pzblockdiag.o ../Matrix/pzbndmat.o ../Matrix/pzespmat.o ../Matrix/pzfmatrix.o \
	../Matrix/pzmatred.o ../Matrix/pzmatrix.o ../Matrix/pzsbndmat.o ../Matrix/pzseqsolver.o \
	../Matrix/pzsespmat.o ../Matrix/pzsfulmat.o ../Matrix/pzshtmat.o ../Matrix/pzskylmat.o \
	../Matrix/pzskylmatpar.o ../Matrix/pzsolve.o ../Matrix/pzstencil.o ../Matrix/pzstepsolver.o \
	../Matrix/pzsysmp.o ../Matrix/pztempmat.o ../Matrix/pzysmp.o ../Matrix/pzmatrixid.o ../Matrix/pzspblockdiagpivot.o \
	../Geom/pzgeopoint.o ../Geom/pzgeoprism.o ../Geom/pzgeopyramid.o ../Geom/pzgeoquad.o \
	../Geom/pzgeotetrahedra.o ../Geom/pzgeotriangle.o ../Geom/TPZGeoCube.o ../Geom/TPZGeoElement.o \
	../Geom/TPZGeoLinear.o ../StrMatrix/pzbdstrmatrix.o ../StrMatrix/pzbstrmatrix.o \
	../StrMatrix/pzfstrmatrix.o ../StrMatrix/pzsbstrmatrix.o ../StrMatrix/pzskylstrmatrix.o \
	../StrMatrix/pzstrmatrix.o ../StrMatrix/TPBSpStructMatrix.o ../StrMatrix/TPZFrontStructMatrix.o \
	../StrMatrix/TPZParFrontStructMatrix.o ../StrMatrix/TPZParSkylineStructMatrix.o \
	../StrMatrix/TPZSpStructMatrix.o ../Frontal/pzvecfrontal.o ../Frontal/tpzeqnarray.o \
	../Frontal/TPZFileEqnStorage.o ../Frontal/TPZFront.o ../Frontal/TPZFrontMatrix.o \
	../Frontal/TPZFrontNonSym.o ../Frontal/TPZFrontSym.o ../Frontal/TPZParFrontMatrix.o \
	../Frontal/TPZStackEqnStorage.o ../Analysis/pzanalysis.o ../Analysis/pzanalysiserror.o \
	../Analysis/pzmganalysis.o ../Analysis/pznonlinanalysis.o ../Analysis/pzonedref.o \
	../Analysis/pzsmanal.o ../Analysis/pzsmfrontalanal.o ../Analysis/pzeuleranalysis.o \
	../Cosys/pzcartsys.o ../Cosys/pzcosys.o ../Cosys/pzcylinsys.o ../Cosys/pzelgq2dcyl.o \
	../Cosys/pzesfersys.o ../Integral/pzquad.o ../Material/eulerdif.o ../Material/pzbndcond.o \
	../Material/pzelasmat.o ../Material/pzmat1dlin.o ../Material/pzmat2dlin.o \
	../Material/pzmaterial.o ../Material/tpzoutofrange.o ../Material/pzmathyperelastic.o \
	../Material/pzmatorthotropic.o ../Material/pzmatplaca2.o ../Material/pzmatplaca.o \
	../Material/pzmattest3d.o ../Material/pzmattest.o ../Material/pzmultplaca.o ../Material/pzplaca.o \
	../Material/pzpoisson3d.o ../Material/TPZEuler.o ../Material/TPZEulerConsLaw.o \
	../Material/TPZLinearConvection.o ../Material/tpzmultcamada.o ../Material/pzbctension.o \
	../Material/TPZMulticamadaOrtho.o ../Material/TPZPlacaOrthotropic.o ../Material/TPZConservationLaw.o \
	../Material/TPZDiffusionConsLaw.o ../Material/pzconslaw.o ../Material/pzartdiff.o \
	../Material/pzeulerconslaw.o ../Material/pzbiharmonic.o ../Material/pzmaterialid.o \
	../Material/tpzoutofrange.o ../Post/pzdxmesh.o ../Post/pzgraphel1d.o ../Post/pzgraphel1dd.o \
	../Post/pzgraphel.o ../Post/pzgraphelq2d.o ../Post/pzgraphelq2dd.o ../Post/pzgraphelq3dd.o \
	../Post/pzgraphmesh.o ../Post/pzgraphnode.o ../Post/pzmvmesh.o ../Post/pztrigraph.o \
	../Post/pztrigraphd.o ../Post/pzv3dmesh.o ../Post/pzvisualmatrix.o ../Multigrid/pzmgsolver.o \
	../Multigrid/pztransfer.o ../Multigrid/pztrnsform.o ../External/pzmetis.o \
	../External/pzrenumbering.o ../Pre/pzgengrid.o ../Pre/pzpargrid.o ../Pre/TPZGMSHReadMesh.o \
	../Pre/TPZExtendGridDimension.o ../Save/pzsave.o ../Save/pzfilebuffer.o ../External/sloan/diamtr.o \
	../External/sloan/gegra.o ../External/sloan/iorele.o ../External/sloan/isorti.o \
	../External/sloan/label.o ../External/sloan/number.o ../External/sloan/profi1.o \
	../External/sloan/rootls.o ../External/sloan/vsrtp1.o





## --| NeoPZ |------------------------------------------------------------------
