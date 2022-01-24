//$Id: pzpostprocanalysis.cpp,v 1.10 2010-11-23 18:58:35 diogo Exp $
#include "TPZLinearAnalysis.h"
#include "pzpostprocanalysis.h"
#include "pzpostprocmat.h"
#include "pzcompelpostproc.h"
#include "pzcmesh.h"
#include "pzgmesh.h"
#include "pzvec.h"
#include "tpzautopointer.h"

#include "pzstring.h"
//#include "pzelastoplasticanalysis.h"
#include "pzcreateapproxspace.h"
#include "pzmultiphysicselement.h" 

#include <map>
#include <set>
#include <stdio.h>
#include "pzlog.h"

#include "pzintel.h"

#include "pzrefpoint.h"
#include "pzgeopoint.h"
#include "pzshapepoint.h"
#include "tpzpoint.h"

#include "pzshapelinear.h"
#include "TPZGeoLinear.h"
#include "TPZRefLinear.h"
#include "tpzline.h"

#include "pzshapetriang.h"
#include "pzreftriangle.h"
#include "pzgeotriangle.h"
#include "tpztriangle.h"

#include "pzrefquad.h"
#include "pzshapequad.h"
#include "pzgeoquad.h"
#include "tpzquadrilateral.h"

#include "pzshapeprism.h"
#include "pzrefprism.h"
#include "pzgeoprism.h"
#include "tpzprism.h"

#include "pzshapetetra.h"
#include "pzreftetrahedra.h"
#include "pzgeotetrahedra.h"
#include "tpztetrahedron.h"

#include "pzshapepiram.h"
#include "pzrefpyram.h"
#include "pzgeopyramid.h"
#include "tpzpyramid.h"

#include "TPZGeoCube.h"
#include "pzshapecube.h"
#include "TPZRefCube.h"
#include "tpzcube.h"
#include "pzelctemp.h"

#ifdef PZ_LOG
static TPZLogger PPAnalysisLogger("pz.analysis.postproc");
#endif

using namespace std;

TPZPostProcAnalysis::TPZPostProcAnalysis() : TPZRegisterClassId(&TPZPostProcAnalysis::ClassId),
fpMainMesh(NULL)
{	
}

TPZPostProcAnalysis::TPZPostProcAnalysis(TPZCompMesh * pRef):TPZRegisterClassId(&TPZPostProcAnalysis::ClassId),
TPZLinearAnalysis(), fpMainMesh(pRef)
{
    
    SetCompMesh(pRef);
    
}

TPZPostProcAnalysis::TPZPostProcAnalysis(const TPZPostProcAnalysis &copy) : TPZRegisterClassId(&TPZPostProcAnalysis::ClassId),
TPZLinearAnalysis(copy), fpMainMesh(0)
{
    
}

TPZPostProcAnalysis &TPZPostProcAnalysis::operator=(const TPZPostProcAnalysis &copy)
{
    SetCompMesh(0);
    return *this;
}

TPZPostProcAnalysis::~TPZPostProcAnalysis()
{
    if (fCompMesh) {
        delete fCompMesh;
    }
    
}

/// Set the computational mesh we are going to post process
void TPZPostProcAnalysis::SetCompMesh(TPZCompMesh *pRef, bool mustOptimizeBandwidth)
{
    // the postprocess mesh already exists, do nothing
    if (fpMainMesh == pRef) {
        return;
    }
    
    if (fCompMesh) {
        delete fCompMesh;
        fCompMesh = 0;
        TPZLinearAnalysis::CleanUp();
    }

    fpMainMesh = pRef;
    
    if (!pRef) {
        return;
    }
    
    TPZCompMesh* pcMainMesh = fpMainMesh;
    
    TPZGeoMesh * pgmesh = pcMainMesh->Reference();

    TPZCompMesh * pcPostProcMesh = new TPZCompMesh(pgmesh);
    
    fCompMesh = pcPostProcMesh;
    
    TPZPostProcAnalysis::SetAllCreateFunctionsPostProc(pcPostProcMesh);
    
}


void TPZPostProcAnalysis::SetPostProcessVariables(TPZVec<int> & matIds, TPZVec<std::string> &varNames)
{

    int nMat, matNumber;
	TPZCompMesh * pcMainMesh = fpMainMesh;

	TPZCompMesh * pcPostProcMesh = (this->Mesh());
    
    if (!pcPostProcMesh) {
        DebugStop();
    }
	
	nMat = matIds.NElements();
	for(int i = 0; i < nMat; i++)
	{
		TPZMaterial * pmat = pcMainMesh->FindMaterial(matIds[i]);
		if(!pmat)
		{
			PZError << "Error at " << __PRETTY_FUNCTION__ << " TPZPostProcAnalysis::SetPostProcessVariables() material Id " << matIds[i] << " not found in original mesh!\n";
			continue;
		}
		
		TPZPostProcMat * pPostProcMat = new TPZPostProcMat(matIds[i]);
		
		pPostProcMat->SetPostProcessVarIndexList(varNames,pmat);
		
		matNumber = pcPostProcMesh->InsertMaterialObject(pPostProcMat);

	}
    
	AutoBuildDisc();
	
    pcPostProcMesh->ExpandSolution();
}

void TPZPostProcAnalysis::AutoBuildDisc() 
{
	TPZAdmChunkVector<TPZGeoEl *> &elvec = Mesh()->Reference()->ElementVec();
	int64_t i, nelem = elvec.NElements();
	int neltocreate = 0;
    
    // build a data structure indicating which geometric elements will be post processed
    fpMainMesh->LoadReferences();
    std::map<TPZGeoEl *,TPZCompEl *> geltocreate;
    TPZCompMesh * pcPostProcMesh = this->Mesh();
    for (i=0; i<nelem; i++) {
        TPZGeoEl * gel = elvec[i];
        if (!gel) {
            continue;
        }
        
        if (gel->HasSubElement()) {
            continue;
        }
        
        TPZMaterial * mat = pcPostProcMesh->FindMaterial(gel->MaterialId());
        if(!mat)
        {
            continue;
        }
        
        if (gel->Reference()) {
            geltocreate[elvec[i]] = gel->Reference();
        }
    }
    Mesh()->Reference()->ResetReference();
    Mesh()->LoadReferences();
    neltocreate = geltocreate.size();
    
	std::set<int> matnotfound;
	int nbl = Mesh()->Block().NBlocks();
	if(neltocreate > nbl) Mesh()->Block().SetNBlocks(neltocreate);
	Mesh()->Block().SetNBlocks(nbl);
    
    std::map<TPZGeoEl *, TPZCompEl *>::iterator it;
    for (it=geltocreate.begin(); it!= geltocreate.end(); it++) 
    {
		TPZGeoEl *gel = it->first;
		if(!gel) continue;
        int matid = gel->MaterialId();
        TPZMaterial * mat = Mesh()->FindMaterial(matid);
        if(!mat)
        {
            matnotfound.insert(matid);
            continue;
        }
        TPZCompEl *cel = Mesh()->CreateCompEl(gel);
        TPZCompElPostProcBase *celpost = dynamic_cast<TPZCompElPostProcBase *>(cel);
        if(!celpost) DebugStop();
        TPZCompEl *celref = it->second;
        int nc = cel->NConnects();
        int ncref = celref->NConnects();
        TPZInterpolationSpace *celspace = dynamic_cast<TPZInterpolationSpace *>(cel);
        TPZInterpolationSpace *celrefspace = dynamic_cast<TPZInterpolationSpace *>(celref);
        int porder;
        if (!celrefspace) {
            TPZMultiphysicsElement *celrefmf = dynamic_cast<TPZMultiphysicsElement *>(celref);
            if (celrefmf){
                celrefspace = dynamic_cast<TPZInterpolationSpace *>(celrefmf->Element(0));
            } else {
                DebugStop();
            }
        }
        celpost->fReferredElement = celrefspace;

        if (celrefspace) {
            porder = celrefspace->GetPreferredOrder();
        } else {
            DebugStop();
        }

        celspace->SetPreferredOrder(porder);
        for (int ic=0; ic<nc; ic++) {
            cel->Connect(ic).SetOrder(porder,cel->ConnectIndex(ic));
            int nshape = celspace->NConnectShapeF(ic,porder);
            cel->Connect(ic).SetNShape(nshape);
        }
        
        TPZIntPoints &intrule = celspace->GetIntegrationRule();
        const TPZIntPoints &intruleref = celref->GetIntegrationRule();
        TPZIntPoints * cloned_rule = intruleref.Clone();
        cel->SetIntegrationRule(cloned_rule);
        
#ifdef PZDEBUG
        if (cel->GetIntegrationRule().NPoints() != intruleref.NPoints()) {
            DebugStop();
        }
#endif
        // this is why the mesh will be discontinuous!!
        gel->ResetReference();
			
	}
    
    // we changed the properties of the connects
    // now synchronize the connect properties with the block sizes
	int64_t nc= Mesh()->NConnects();
    for (int64_t ic=0; ic<nc; ic++) {
        TPZConnect &c = Mesh()->ConnectVec()[ic];
        int blsize = c.NShape()*c.NState();
        int64_t seqnum = c.SequenceNumber();
        Mesh()->Block().Set(seqnum, blsize);
    }
	Mesh()->InitializeBlock();
#ifdef PZ_LOG
    if(PPAnalysisLogger.isDebugEnabled())
    {
        std::stringstream sout;
        Mesh()->Print(sout);
        LOGPZ_DEBUG(PPAnalysisLogger, sout.str())
    }
#endif
    
#ifdef PZDEBUG
	if(matnotfound.size())
	{
		std::cout << "Post-processing mesh was created without these materials: ";
		std::set<int>::iterator it;
		for(it = matnotfound.begin(); it!= matnotfound.end(); it++)
		{
			std::cout << *it << " ";
		}
		std::cout << std::endl;
        DebugStop();
	}
#endif
	
}

void TPZPostProcAnalysis::Assemble()
{
   PZError << "Error at " << __PRETTY_FUNCTION__ << " TPZPostProcAnalysis::Assemble() should never be called\n";
}

void TPZPostProcAnalysis::Solve(){
   PZError << "Error at " << __PRETTY_FUNCTION__ << " TPZPostProcAnalysis::Solve() should never be called\n";
}

void TPZPostProcAnalysis::TransferSolution()
{

    // this is where we compute the projection of the post processed variables
    TPZLinearAnalysis::AssembleResidual();
    fSolution = Rhs();
    TPZLinearAnalysis::LoadSolution();
    
    TPZCompMesh *compmeshPostProcess = (Mesh());
    if (!compmeshPostProcess) {
        DebugStop();
    }
    // fpMainMesh is the mesh with the actual finite element approximation, but probably stored at
    // integration points
    TPZCompMesh *solmesh = fpMainMesh;
    fpMainMesh->Reference()->ResetReference();
    fpMainMesh->LoadReferences();
    //In case the post processing computed element solutions
    // copy the values from the post processing mesh to the finite element mesh
    TPZFMatrix<STATE> &comprefElSol = compmeshPostProcess->ElementSolution();
    // solmesh if the finite element simulation mesh
    const TPZFMatrix<STATE> &solmeshElSol = solmesh->ElementSolution();
    int64_t numelsol = solmesh->ElementSolution().Cols();
    int64_t nelem = compmeshPostProcess->NElements();
    compmeshPostProcess->ElementSolution().Redim(nelem, numelsol);
    if (numelsol) 
    {
        for (int64_t el=0; el<nelem; el++) {
            TPZCompEl *celpost = compmeshPostProcess->Element(el);
            TPZGeoEl *gel = celpost->Reference();
            // we dont acount for condensed elements submeshes etc
            if(!gel) DebugStop();
            TPZCompEl *cel = gel->Reference();
            if (!cel) {
                DebugStop();
            }
            int64_t index = cel->Index();
            // we copy from the simulation mesh to the post processing mesh
            for (int64_t isol=0; isol<numelsol; isol++) {
                comprefElSol(el,isol) = solmeshElSol.Get(index,isol);
            }
        }
    }
}


void TPZPostProcAnalysis::SetAllCreateFunctionsPostProc(TPZCompMesh *cmesh)
{
    
    TPZManVector<TCreateFunction,10> functions(8);
    
    functions[EPoint] = &TPZPostProcAnalysis::CreatePointEl;
    functions[EOned] = TPZPostProcAnalysis::CreateLinearEl;
    functions[EQuadrilateral] = TPZPostProcAnalysis::CreateQuadEl;
    functions[ETriangle] = TPZPostProcAnalysis::CreateTriangleEl;
    functions[EPrisma] = TPZPostProcAnalysis::CreatePrismEl;
    functions[ETetraedro] = TPZPostProcAnalysis::CreateTetraEl;
    functions[EPiramide] = TPZPostProcAnalysis::CreatePyramEl;
    functions[ECube] = TPZPostProcAnalysis::CreateCubeEl;
    cmesh->ApproxSpace().SetCreateFunctions(functions);
}


#include "TPZCompElH1.h"

using namespace pzshape;

template class TPZCompElPostProc< TPZCompElH1<TPZShapePoint> >;
template class TPZCompElPostProc< TPZCompElH1<TPZShapeLinear> >;
template class TPZCompElPostProc< TPZCompElH1<TPZShapeQuad> >;
template class TPZCompElPostProc< TPZCompElH1<TPZShapeTriang> >;
template class TPZCompElPostProc< TPZCompElH1<TPZShapeCube> >;
template class TPZCompElPostProc< TPZCompElH1<TPZShapePrism> >;
template class TPZCompElPostProc< TPZCompElH1<TPZShapePiram> >;
template class TPZCompElPostProc< TPZCompElH1<TPZShapeTetra> >;
template class TPZCompElPostProc< TPZCompElDisc >;

template class TPZRestoreClass<TPZCompElPostProc< TPZCompElH1<TPZShapePoint> >>;
template class TPZRestoreClass<TPZCompElPostProc< TPZCompElH1<TPZShapeLinear> >>;
template class TPZRestoreClass<TPZCompElPostProc< TPZCompElH1<TPZShapeQuad> >>;
template class TPZRestoreClass<TPZCompElPostProc< TPZCompElH1<TPZShapeTriang> >>;
template class TPZRestoreClass<TPZCompElPostProc< TPZCompElH1<TPZShapeCube> >>;
template class TPZRestoreClass<TPZCompElPostProc< TPZCompElH1<TPZShapePrism> >>;
template class TPZRestoreClass<TPZCompElPostProc< TPZCompElH1<TPZShapePiram> >>;
template class TPZRestoreClass<TPZCompElPostProc< TPZCompElH1<TPZShapeTetra> >>;
template class TPZRestoreClass<TPZCompElPostProc< TPZCompElDisc >>;

TPZCompEl *TPZPostProcAnalysis::CreatePointEl(TPZGeoEl *gel,TPZCompMesh &mesh) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElPostProc< TPZCompElH1<TPZShapePoint> >(mesh,gel);
	return NULL;
}
TPZCompEl *TPZPostProcAnalysis::CreateLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElPostProc<TPZCompElH1<TPZShapeLinear> >(mesh,gel);
	return NULL;
}
TPZCompEl *TPZPostProcAnalysis::CreateQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElPostProc<TPZCompElH1<TPZShapeQuad> >(mesh,gel);
	return NULL;
}
TPZCompEl *TPZPostProcAnalysis::CreateTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElPostProc<TPZCompElH1<TPZShapeTriang> >(mesh,gel);
	return NULL;
}
TPZCompEl *TPZPostProcAnalysis::CreateCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElPostProc<TPZCompElH1<TPZShapeCube> >(mesh,gel);
	return NULL;
}
TPZCompEl *TPZPostProcAnalysis::CreatePrismEl(TPZGeoEl *gel,TPZCompMesh &mesh) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElPostProc< TPZCompElH1<TPZShapePrism> >(mesh,gel);
	return NULL;
}
TPZCompEl *TPZPostProcAnalysis::CreatePyramEl(TPZGeoEl *gel,TPZCompMesh &mesh) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElPostProc<TPZCompElH1<TPZShapePiram> >(mesh,gel);
	return NULL;
}
TPZCompEl *TPZPostProcAnalysis::CreateTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElPostProc<TPZCompElH1<TPZShapeTetra> >(mesh,gel);
	return NULL;
}


TPZCompEl * TPZPostProcAnalysis::CreatePostProcDisc(TPZGeoEl *gel, TPZCompMesh &mesh)
{
	return new TPZCompElPostProc< TPZCompElDisc > (mesh,gel);
}

/** @brief Returns the unique identifier for reading/writing objects to streams */
int TPZPostProcAnalysis::ClassId() const{
    return Hash("TPZPostProcAnalysis") ^ TPZLinearAnalysis::ClassId() << 1;
}
/** @brief Save the element data to a stream */
void TPZPostProcAnalysis::Write(TPZStream &buf, int withclassid) const
{
    TPZLinearAnalysis::Write(buf, withclassid);
    TPZPersistenceManager::WritePointer(fpMainMesh, &buf);
}

/** @brief Read the element data from a stream */
void TPZPostProcAnalysis::Read(TPZStream &buf, void *context)
{
    TPZLinearAnalysis::Read(buf, context);
    fpMainMesh = dynamic_cast<TPZCompMesh*>(TPZPersistenceManager::GetInstance(&buf));
}

