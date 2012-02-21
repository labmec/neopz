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

#include <map>
#include <set>
#include <stdio.h>
#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr PPAnalysisLogger(Logger::getLogger("analysis.postproc"));
#endif

using namespace std;

TPZPostProcAnalysis::TPZPostProcAnalysis() : fpMainAnalysis(NULL)
{	
}

TPZPostProcAnalysis::TPZPostProcAnalysis(TPZAnalysis * pRef):TPZAnalysis(), fpMainAnalysis(pRef)
{
	
	if(!fpMainAnalysis)
	{
		PZError << "Error at " << __PRETTY_FUNCTION__ << " TPZPostProcAnalysis::TPZPostProcAnalysis() Invalid analysis to post process!\n";
		return;
	}

	TPZCompMesh* pcMainMesh = fpMainAnalysis->Mesh();
	//pcMainMesh->SetAllCreateFunctionsContinuous();
	if(!pcMainMesh)
	{
		PZError << "Error at " << __PRETTY_FUNCTION__ << " TPZPostProcAnalysis::TPZPostProcAnalysis() Invalid fCompMesh in analysis to post process!\n";
		return;
	}
	
	TPZGeoMesh * pgmesh = pcMainMesh->Reference();
	
	if(!pgmesh)
	{
		PZError << "Error at " << __PRETTY_FUNCTION__ << " TPZPostProcAnalysis::TPZPostProcAnalysis() Invalid GeoMesh in analysis to post process!\n";
		return;
	}
	
	

	TPZCompMeshReferred * pcPostProcMesh = new TPZCompMeshReferred(pgmesh);
	SetAllCreateFunctionsPostProc(pcPostProcMesh);
	
	
	fCompMesh = pcPostProcMesh;
}

TPZPostProcAnalysis::~TPZPostProcAnalysis()
{
}

void TPZPostProcAnalysis::SetPostProcessVariables(TPZVec<int> & matIds, TPZVec<std::string> &varNames)
{
	int i, j, nMat, matNumber;

	if(!fpMainAnalysis)
	{
		PZError << "Error at " << __PRETTY_FUNCTION__ << " TPZPostProcAnalysis::SetPostProcessVariables() unavailable fpRef!\n";
		return;	
	}

	TPZCompMesh * pcMainMesh = fpMainAnalysis->Mesh();
	if(!pcMainMesh)
	{
		PZError << "Error at " << __PRETTY_FUNCTION__ << " TPZPostProcAnalysis::SetPostProcessVariables() unavailable fpRef->Mesh()!\n";
		return;
	}
	
	TPZGeoMesh * pgmesh = pcMainMesh->Reference();
	if(!pgmesh)
	{
		PZError << "Error at " << __PRETTY_FUNCTION__ << " TPZPostProcAnalysis::TPZPostProcAnalysis() Invalid GeoMesh in analysis to post process!\n";
		return;
	}
	
	TPZCompMeshReferred * pcPostProcMesh = dynamic_cast<TPZCompMeshReferred *>(this->Mesh());
	if(!pcPostProcMesh)
	{
		PZError << "Error at " << __PRETTY_FUNCTION__ << " TPZPostProcAnalysis::SetPostProcessVariables() unavailable this->Mesh()!\n";
		return;
	}
	
	TPZStack<int> avlMatIds;
	int nel = pgmesh->NElements();	
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
	
	nMat = matIds.NElements();
	for(i = 0; i < nMat; i++)
	{
		TPZAutoPointer<TPZMaterial> pmat = pcMainMesh->FindMaterial(matIds[i]);
		if(!pmat)
		{
			PZError << "Error at " << __PRETTY_FUNCTION__ << " TPZPostProcAnalysis::SetPostProcessVariables() material Id " << matIds[i] << " not found in original mesh!\n";
			continue;
		}
		
		TPZPostProcMat * pPostProcMat = new TPZPostProcMat(matIds[i]);
		
		pPostProcMat->SetPostProcessVarIndexList(varNames,pmat.operator->());
		
		matNumber = pcPostProcMesh->InsertMaterialObject(TPZAutoPointer<TPZMaterial>(pPostProcMat));
	}
	
//	pcPostProcMesh->AutoBuild();
	//pcPostProcMesh->AutoBuildDisc();
	AutoBuildDisc();
	
	pcPostProcMesh->LoadReferred(pcMainMesh);
}

void TPZPostProcAnalysis::AutoBuildDisc() 
{
	TPZAdmChunkVector<TPZGeoEl *> &elvec = Mesh()->Reference()->ElementVec();
	int i, nelem = elvec.NElements();
	int neltocreate = 0;
	int index;
	for(i=0; i<nelem; i++) {
		TPZGeoEl *gel = elvec[i];
		if(!gel) continue;
		if(!gel->HasSubElement()) {
			neltocreate++;
		}
	}
	std::set<int> matnotfound;
	int nbl = Mesh()->Block().NBlocks();
	if(neltocreate > nbl) Mesh()->Block().SetNBlocks(neltocreate);
	Mesh()->Block().SetNBlocks(nbl);
	for(i=0; i<nelem; i++) {
		TPZGeoEl *gel = elvec[i];
		if(!gel) continue;
		if(!gel->HasSubElement()) {
			int matid = gel->MaterialId();
			TPZAutoPointer<TPZMaterial> mat = Mesh()->FindMaterial(matid);
			if(!mat)
			{
				matnotfound.insert(matid);
				continue;
			}
			int printing = 1;
			if (printing) {
				gel->Print(cout);
			}
			
			//if(!gel->Reference() && gel->NumInterfaces() == 0)
			Mesh()->CreateCompEl(gel,index);
			gel->ResetReference();
			
		}
	}
	
	Mesh()->InitializeBlock();
	if(matnotfound.size())
	{
		std::cout << "Malha post proc was not created properly because of these materials ";
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

#include "pzcreateapproxspace.h"



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
using namespace pzshape;

TPZCompEl *TPZPostProcAnalysis::CreatePointEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElPostProc< TPZIntelGen<TPZShapePoint> >(mesh,gel,index);
	return NULL;
}
TPZCompEl *TPZPostProcAnalysis::CreateLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElPostProc<TPZIntelGen<TPZShapeLinear> >(mesh,gel,index);
	return NULL;
}
TPZCompEl *TPZPostProcAnalysis::CreateQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElPostProc<TPZIntelGen<TPZShapeQuad> >(mesh,gel,index);
	return NULL;
}
TPZCompEl *TPZPostProcAnalysis::CreateTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElPostProc<TPZIntelGen<TPZShapeTriang> >(mesh,gel,index);
	return NULL;
}
TPZCompEl *TPZPostProcAnalysis::CreateCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElPostProc<TPZIntelGen<TPZShapeCube> >(mesh,gel,index);
	return NULL;
}
TPZCompEl *TPZPostProcAnalysis::CreatePrismEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElPostProc< TPZIntelGen<TPZShapePrism> >(mesh,gel,index);
	return NULL;
}
TPZCompEl *TPZPostProcAnalysis::CreatePyramEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElPostProc<TPZIntelGen<TPZShapePiram> >(mesh,gel,index);
	return NULL;
}
TPZCompEl *TPZPostProcAnalysis::CreateTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	if(!gel->Reference() && gel->NumInterfaces() == 0)
		return new TPZCompElPostProc<TPZIntelGen<TPZShapeTetra> >(mesh,gel,index);
	return NULL;
}


TPZCompEl * TPZPostProcAnalysis::CreatePostProcDisc(TPZGeoEl *gel, TPZCompMesh &mesh, int &index)
{
	return new TPZCompElPostProc< TPZCompElDisc > (mesh,gel,index);
}
