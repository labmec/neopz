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

libpz_a_LIBADD = ../Post/tpzgraphelt2dmapped.o \
	../Post/tpzgraphelpyramidmapped.o ../Post/tpzgraphelprismmapped.o ../Analysis/TPZNLMultGridAnalysis.o \
	../Analysis/pzanalysis.o ../Analysis/pzanalysiserror.o ../Analysis/pzeuleranalysis.o \
	../Analysis/pzmganalysis.o ../Analysis/pznonlinanalysis.o ../Analysis/pzsmanal.o \
	../Analysis/pzsmfrontalanal.o ../Analysis/pztransientanalysis.o ../Common/pzreal.o \
	../Cosys/pzcartsys.o ../Cosys/pzcosys.o ../Cosys/pzcylinsys.o ../Cosys/pzesfersys.o \
	../External/TPZBoostGraph.o ../External/pzmetis.o ../External/pzrenumbering.o ../External/pzsloan.o \
	../External/sloan/diamtr.o ../External/sloan/gegra.o ../External/sloan/iorele.o \
	../External/sloan/isorti.o ../External/sloan/label.o ../External/sloan/number.o \
	../External/sloan/profi1.o ../External/sloan/rootls.o ../External/sloan/vsrtp1.o \
	../External/tpznodesetcompute.o ../Frontal/TPZFileEqnStorage.o ../Frontal/TPZFront.o \
	../Frontal/TPZFrontMatrix.o ../Frontal/TPZFrontNonSym.o ../Frontal/TPZFrontSym.o \
	../Frontal/TPZParFrontMatrix.o ../Frontal/TPZStackEqnStorage.o ../Frontal/pzvecfrontal.o \
	../Frontal/tpzeqnarray.o ../Geom/TPZGeoCube.o ../Geom/TPZGeoElement.o ../Geom/TPZGeoLinear.o \
	../Geom/pzgeopoint.o ../Geom/pzgeoprism.o ../Geom/pzgeopyramid.o ../Geom/pzgeoquad.o \
	../Geom/pzgeotetrahedra.o ../Geom/pzgeotriangle.o ../Geom/pznoderep.o ../Integral/pzquad.o \
	../Material/TPZConservationLaw.o ../Material/TPZDiffusionConsLaw.o ../Material/TPZEuler.o \
	../Material/TPZEulerConsLaw.o ../Material/TPZLinearConvection.o ../Material/TPZMulticamadaOrtho.o \
	../Material/TPZPlacaOrthotropic.o ../Material/eulerdif.o ../Material/pzartdiff.o \
	../Material/pzbctension.o ../Material/pzbiharmonic.o ../Material/pzbndcond.o \
	../Material/pzconslaw.o ../Material/pzcoupledtransportdarcy.o \
	../Material/pzcoupledtransportdarcyBC.o ../Material/pzdiscgal.o ../Material/pzelasmat.o ../Material/pzelast3d.o \
	../Material/pzeulerconslaw.o ../Material/pzincnskeps.o ../Material/pzl2projection.o \
	../Material/pzmat1dlin.o ../Material/pzmat2dlin.o ../Material/pzmaterial.o \
	../Material/pzmaterialdata.o ../Material/pzmaterialid.o ../Material/pzmathyperelastic.o \
	../Material/pzmatorthotropic.o ../Material/pzmatplaca2.o \
	../Material/pzmattest.o ../Material/pzmattest3d.o ../Material/pzmultplaca.o \
	../Material/pznlmat1d.o ../Material/pznlmat1drotatedengstrain.o \
	../Material/pznonlinbiharmonic.o ../Material/pznonlinearpoisson3d.o ../Material/pzplaca.o \
	../Material/pzpoisson3d.o ../Material/pzpoisson3dreferred.o ../Material/pzthermicelast3d.o \
	../Material/pztransientmat.o ../Material/tpzmultcamada.o ../Material/tpzoutofrange.o \
	../Matrix/pzblock.o ../Matrix/pzblockdiag.o ../Matrix/pzbndmat.o ../Matrix/pzespmat.o \
	../Matrix/pzfmatrix.o ../Matrix/pzlink.o ../Matrix/pzmatred.o ../Matrix/pzmatrix.o \
	../Matrix/pzmatrixid.o ../Matrix/pzsbndmat.o ../Matrix/pzseqsolver.o ../Matrix/pzsespmat.o \
	../Matrix/pzsfulmat.o ../Matrix/pzshtmat.o ../Matrix/pzskylmat.o ../Matrix/pzskylmatpar.o \
	../Matrix/pzsolve.o ../Matrix/pzspblockdiagpivot.o ../Matrix/pzstencil.o \
	../Matrix/pzstepsolver.o ../Matrix/pzsysmp.o ../Matrix/pzysmp.o \
	../Matrix/tpzsparseblockdiagonal.o ../Mesh/TPZAgglomerateEl.o ../Mesh/TPZCompElDisc.o \
	../Mesh/TPZInterfaceEl.o ../Mesh/pzcheckgeom.o ../Mesh/pzcheckmesh.o ../Mesh/pzcheckrestraint.o \
	../Mesh/pzcmesh.o ../Mesh/pzcompel.o ../Mesh/pzconnect.o ../Mesh/pzelctemp.o \
	../Mesh/pzelmat.o ../Mesh/pzflowcmesh.o ../Mesh/pzgeoel.o ../Mesh/pzgeoelbc.o \
	../Mesh/pzgeoelrefless.o ../Mesh/pzgeoelside.o ../Mesh/pzgmesh.o ../Mesh/pzgnode.o \
	../Mesh/pzintel.o ../Mesh/pzinterpolationspace.o ../Mesh/pzmeshid.o \
	../Mesh/pzreferredcompel.o ../Mesh/pzstackmesh.o ../Mesh/pzsubcmesh.o ../Mesh/pzvecmesh.o \
	../Mesh/tpzcompmeshreferred.o ../Mesh/tpzgeoelrefpattern.o ../Multigrid/pzmgsolver.o \
	../Multigrid/pztransfer.o ../Multigrid/pztrnsform.o ../Post/pzdxmesh.o ../Post/pzgraphel.o \
	../Post/pzgraphel1d.o ../Post/pzgraphel1dd.o ../Post/pzgraphelq2d.o ../Post/pzgraphelq2dd.o \
	../Post/pzgraphelq3dd.o ../Post/pzgraphmesh.o ../Post/pzgraphnode.o ../Post/pzmvmesh.o \
	../Post/pztrigraph.o ../Post/pztrigraphd.o ../Post/pzv3dmesh.o ../Post/pzvisualmatrix.o \
	../Post/tpzgraphelt3d.o ../Pre/TPZExtendGridDimension.o ../Pre/TPZGMSHReadMesh.o \
	../Pre/pzgengrid.o ../Pre/pzpargrid.o ../Pre/pzreadmesh.o ../Pre/pzreadmeshhr.o \
	../Pre/pzreadtetgen.o ../Refine/TPZRefCube.o ../Refine/TPZRefLinear.o \
	../Refine/TPZRefPattern.o ../Refine/pzrefpoint.o ../Refine/pzrefprism.o ../Refine/pzrefpyram.o \
	../Refine/pzrefquad.o ../Refine/pzreftetrahedra.o ../Refine/pzreftriangle.o \
	../Save/pzfilebuffer.o ../Save/pzsave.o ../Shape/TPZShapeDisc.o ../Shape/pzshapecube.o \
	../Shape/pzshapelinear.o ../Shape/pzshapepiram.o ../Shape/pzshapeprism.o ../Shape/pzshapequad.o \
	../Shape/pzshapetetra.o ../Shape/pzshapetriang.o ../StrMatrix/TPBSpStructMatrix.o \
	../StrMatrix/TPZFrontStructMatrix.o ../StrMatrix/TPZParFrontStructMatrix.o \
	../StrMatrix/TPZParSkylineStructMatrix.o ../StrMatrix/TPZSpStructMatrix.o ../StrMatrix/pzbdstrmatrix.o \
	../StrMatrix/pzbstrmatrix.o ../StrMatrix/pzfstrmatrix.o ../StrMatrix/pzsbstrmatrix.o \
	../StrMatrix/pzskylstrmatrix.o ../StrMatrix/pzstrmatrix.o ../Topology/tpzcube.o ../Topology/tpzline.o \
	../Topology/tpzpoint.o ../Topology/tpzprism.o ../Topology/tpzpyramid.o \
	../Topology/tpzquadrilateral.o ../Topology/tpztetrahedron.o ../Topology/tpztriangle.o \
	../Util/TPZTimer.o ../Util/pzfunction.o ../Util/pzlog.o ../Util/pzmap.o ../Util/pzstring.o \
	../Util/tpzautopointer.o ../Util/tpzpermutation.o

## --| NeoPZ |------------------------------------------------------------------
