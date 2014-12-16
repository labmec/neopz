//$Id: pzpostprocanalysis.cpp,v 1.10 2010-11-23 18:58:35 diogo Exp $
#include "pzanalysis.h"
#include "pzpostprocanalysis.h"
#include "pzpostprocmat.h"
#include "pzcompelpostproc.h"
#include "pzcmesh.h"
#include "pzgmesh.h"
#include "pzvec.h"
#include "tpzautopointer.h"
#include "tpzcompmeshreferred.h"
#include "pzstring.h"
//#include "pzelastoplasticanalysis.h"
#include "pzcreateapproxspace.h"
#include "pzmeshid.h"

#include <map>
#include <set>
#include <stdio.h>
#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr PPAnalysisLogger(Logger::getLogger("pz.analysis.postproc"));
#endif

using namespace std;

TPZPostProcAnalysis::TPZPostProcAnalysis() : fpMainMesh(NULL)
{	
}

TPZPostProcAnalysis::TPZPostProcAnalysis(TPZCompMesh * pRef):TPZAnalysis(), fpMainMesh(pRef)
{
    
    SetCompMesh(pRef);
    
}

TPZPostProcAnalysis::TPZPostProcAnalysis(const TPZPostProcAnalysis &copy) : TPZAnalysis(), fpMainMesh(0)
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
void TPZPostProcAnalysis::SetCompMesh(TPZCompMesh *pRef)
{
    // the postprocess mesh already exists, do nothing
    if (fpMainMesh == pRef) {
        return;
    }
    
    if (fCompMesh) {
//        std::cout << "PostProcAnalysis deleting the mesh " << (void *) fCompMesh << std::endl;
        delete fCompMesh;
        fCompMesh = 0;
        TPZAnalysis::CleanUp();
    }

    fpMainMesh = pRef;
    
    if (!pRef) {
        return;
    }
    
    TPZCompMesh* pcMainMesh = fpMainMesh;
    
    TPZGeoMesh * pgmesh = pcMainMesh->Reference();
    
    // TPZPostProcAnalysis::SetAllCreateFunctionsPostProc();
    
    
    TPZCompMeshReferred * pcPostProcMesh = new TPZCompMeshReferred(pgmesh);
    
    fCompMesh = pcPostProcMesh;
    
    TPZPostProcAnalysis::SetAllCreateFunctionsPostProc(pcPostProcMesh);
    
}


void TPZPostProcAnalysis::SetPostProcessVariables(TPZVec<int> & matIds, TPZVec<std::string> &varNames)
{
	//int j;
    int nMat, matNumber;


	TPZCompMesh * pcMainMesh = fpMainMesh;
	
	//TPZGeoMesh * pgmesh = pcMainMesh->Reference();

	TPZCompMeshReferred * pcPostProcMesh = dynamic_cast<TPZCompMeshReferred *>(this->Mesh());
    
    if (!pcPostProcMesh) {
        DebugStop();
    }
	
    if (pcPostProcMesh->ReferredMesh() == pcMainMesh) {
        return;
    }

    /*
	TPZStack<int> avlMatIds;
	long nel = pgmesh->NElements(), i;	
	for(i = 0; i < nel; i++)
	{
		int matId = pgmesh->ElementVec()[i]->MaterialId();
		int isMatPostProc = 0;
		int isMatAvl = 0;
		j = 0;
		nMat = matIds.NElements();
		while(j < nMat && !isMatPostProc)
		{
			if(matId == matIds[j])isMatPostProc = 1;
			j++;
		}
		
		if(!isMatPostProc)
		{
			nMat = avlMatIds.NElements();
			j = 0;
			while(j < nMat && !isMatAvl)
			{
				if(matId == avlMatIds[j])isMatAvl = 1;
				j++;
			}
			
			if(!isMatAvl)
			{
				avlMatIds.Push(matId);
			}
		}
	}
	*/
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
	
	//pcPostProcMesh->AutoBuild();
	//pcPostProcMesh->AutoBuildDisc();
    
	AutoBuildDisc();
	
	pcPostProcMesh->LoadReferred(pcMainMesh);
}

void TPZPostProcAnalysis::AutoBuildDisc() 
{
	TPZAdmChunkVector<TPZGeoEl *> &elvec = Mesh()->Reference()->ElementVec();
	long i, nelem = elvec.NElements();
	int neltocreate = 0;
	long index;
    // build a data structure indicating which geometric elements will be post processed
    fpMainMesh->LoadReferences();
    std::map<TPZGeoEl *,TPZCompEl *> geltocreate;
    for (i=0; i<nelem; i++) {
        if (!elvec[i]) {
            continue;
        }
        if (elvec[i]->Reference()) {
            geltocreate[elvec[i]] = elvec[i]->Reference();
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
        int printing = 0;
        if (printing) {
            gel->Print(cout);
        }
			
        Mesh()->CreateCompEl(gel,index);
        TPZCompEl *cel = Mesh()->ElementVec()[index];
        TPZCompEl *celref = it->second;
        int nc = cel->NConnects();
        int ncref = celref->NConnects();
        if (nc != ncref) {
            DebugStop();
        }
        TPZInterpolationSpace *celspace = dynamic_cast<TPZInterpolationSpace *>(cel);
        TPZInterpolationSpace *celrefspace = dynamic_cast<TPZInterpolationSpace *>(celref);
        int porder = celrefspace->GetPreferredOrder();
//        if (porder != 2) {
//            std::cout << "I should stop porder = " << porder << std::endl;
//        }
        celspace->SetPreferredOrder(porder);
        for (int ic=0; ic<nc; ic++) {
            int conorder = celref->Connect(ic).Order();
            cel->Connect(ic).SetOrder(conorder);
            int nshape = celspace->NConnectShapeF(ic);
            cel->Connect(ic).SetNShape(nshape);
        }
        TPZIntPoints &intrule = celspace->GetIntegrationRule();
        TPZVec<int> intorder(gel->Dimension(),0);
        const TPZIntPoints &intruleref = celrefspace->GetIntegrationRule();
        intruleref.GetOrder(intorder);
        intrule.SetOrder(intorder);
#ifdef DEBUG
        if (intrule.NPoints() != intruleref.NPoints()) {
            DebugStop();
        }
#endif
        gel->ResetReference();
			
	}
    
    // we changed the properties of the connects
    // now synchronize the connect properties with the block sizes
	long nc= Mesh()->NConnects();
    for (long ic=0; ic<nc; ic++) {
        TPZConnect &c = Mesh()->ConnectVec()[ic];
        int blsize = c.NShape()*c.NState();
        long seqnum = c.SequenceNumber();
        Mesh()->Block().Set(seqnum, blsize);
    }
	Mesh()->InitializeBlock();
#ifdef LOG4CXX
    if(PPAnalysisLogger->isDebugEnabled())
    {
        std::stringstream sout;
        Mesh()->Print(sout);
        LOGPZ_DEBUG(PPAnalysisLogger, sout.str())
    }
#endif
	if(matnotfound.size())
	{
		std::cout << "Malha post proc was created without these materials ";
		std::set<int>::iterator it;
		for(it = matnotfound.begin(); it!= matnotfound.end(); it++)
		{
			std::cout << *it << " ";
		}
		std::cout << std::endl;
	}
	
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

    
    TPZAnalysis::AssembleResidual();
    fSolution = Rhs();
    TPZAnalysis::LoadSolution();
    
    TPZCompMeshReferred *compref = dynamic_cast<TPZCompMeshReferred *>(Mesh());
    if (!compref) {
        DebugStop();
    }
    TPZCompMesh *solmesh = fpMainMesh;
    long numelsol = solmesh->ElementSolution().Cols();
    long nelem = compref->NElements();
    compref->ElementSolution().Redim(nelem, numelsol);
    if (numelsol) 
    {
        for (long el=0; el<nelem; el++) {
            TPZCompEl *cel = compref->ReferredEl(el);
            if (!cel) {
                continue;
            }
            long index = cel->Index();
            for (long isol=0; isol<numelsol; isol++) {
                compref->ElementSolution()(el,isol) = solmesh->ElementSolution()(index,isol);
            }
        }
    }
}


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


/*void TPZPostProcAnalysis::SetAllCreateFunctionsPostProc()
{
	pzgeom::TPZGeoPoint::fp = TPZPostProcAnalysis::CreatePointEl;
	pzgeom::TPZGeoLinear::fp = TPZPostProcAnalysis::CreateLinearEl;
	pzgeom::TPZGeoQuad::fp = TPZPostProcAnalysis::CreateQuadEl;
	pzgeom::TPZGeoTriangle::fp = TPZPostProcAnalysis::CreateTriangleEl;
	pzgeom::TPZGeoPrism::fp = TPZPostProcAnalysis::CreatePrismEl;
	pzgeom::TPZGeoTetrahedra::fp = TPZPostProcAnalysis::CreateTetraEl;
	pzgeom::TPZGeoPyramid::fp = TPZPostProcAnalysis::CreatePyramEl;
	pzgeom::TPZGeoCube::fp = TPZPostProcAnalysis::CreateCubeEl;


}*/


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
  /*
    functions[EPoint] = TPZCompElDisc::CreateDisc;
    functions[EOned] = TPZCompElDisc::CreateDisc;
    functions[ETriangle] = TPZCompElDisc::CreateDisc;
    functions[EQuadrilateral] = TPZCompElDisc::CreateDisc;
    functions[ETetraedro] = TPZCompElDisc::CreateDisc;
    functions[EPiramide] = TPZCompElDisc::CreateDisc;
    functions[EPrisma] = TPZCompElDisc::CreateDisc;
    functions[ECube] = TPZCompElDisc::CreateDisc;
    cmesh->ApproxSpace().SetCreateFunctions(functions);
   */ 
}



using namespace pzshape;

TPZCompEl *TPZPostProcAnalysis::CreatePointEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElPostProc< TPZIntelGen<TPZShapePoint> >(mesh,gel,index);
	return NULL;
}
TPZCompEl *TPZPostProcAnalysis::CreateLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElPostProc<TPZIntelGen<TPZShapeLinear> >(mesh,gel,index);
	return NULL;
}
TPZCompEl *TPZPostProcAnalysis::CreateQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElPostProc<TPZIntelGen<TPZShapeQuad> >(mesh,gel,index);
	return NULL;
}
TPZCompEl *TPZPostProcAnalysis::CreateTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElPostProc<TPZIntelGen<TPZShapeTriang> >(mesh,gel,index);
	return NULL;
}
TPZCompEl *TPZPostProcAnalysis::CreateCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElPostProc<TPZIntelGen<TPZShapeCube> >(mesh,gel,index);
	return NULL;
}
TPZCompEl *TPZPostProcAnalysis::CreatePrismEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElPostProc< TPZIntelGen<TPZShapePrism> >(mesh,gel,index);
	return NULL;
}
TPZCompEl *TPZPostProcAnalysis::CreatePyramEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElPostProc<TPZIntelGen<TPZShapePiram> >(mesh,gel,index);
	return NULL;
}
TPZCompEl *TPZPostProcAnalysis::CreateTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElPostProc<TPZIntelGen<TPZShapeTetra> >(mesh,gel,index);
	return NULL;
}


TPZCompEl * TPZPostProcAnalysis::CreatePostProcDisc(TPZGeoEl *gel, TPZCompMesh &mesh, long &index)
{
	return new TPZCompElPostProc< TPZCompElDisc > (mesh,gel,index);
}

/** @brief Returns the unique identifier for reading/writing objects to streams */
int TPZPostProcAnalysis::ClassId() const
{
	cout << "\nFIX ME: TPZPostProcAnalysis::ClassId()" << endl;
	DebugStop();
    return -1;
}
/** @brief Save the element data to a stream */
void TPZPostProcAnalysis::Write(TPZStream &buf, int withclassid)
{
    DebugStop();
}

/** @brief Read the element data from a stream */
void TPZPostProcAnalysis::Read(TPZStream &buf, void *context)
{
    DebugStop();
}

