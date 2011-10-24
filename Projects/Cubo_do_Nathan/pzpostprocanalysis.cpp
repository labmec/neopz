//$Id: pzpostprocanalysis.cpp,v 1.6 2010-06-11 22:13:02 diogo Exp $
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


TPZPostProcAnalysis::TPZPostProcAnalysis() : fpMainAnalysis(NULL){
	
}

TPZPostProcAnalysis::TPZPostProcAnalysis(TPZAnalysis * pRef):TPZAnalysis(), fpMainAnalysis(pRef){
	
	if(!fpMainAnalysis)
	{
		PZError << "Error at " << __PRETTY_FUNCTION__ << " TPZPostProcAnalysis::TPZPostProcAnalysis() Invalid analysis to post process!\n";
		return;
	}
	
	TPZCompMesh* pcMainMesh = fpMainAnalysis->Mesh();
	
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
	
	TPZPostProcAnalysis::SetAllCreateFunctionsPostProc();
	
	TPZCompMeshReferred * pcPostProcMesh = new TPZCompMeshReferred(pgmesh);
	
	//pcPostProcMesh->LoadReferred(pcMainMesh);
	
	fCompMesh = pcPostProcMesh;
}

TPZPostProcAnalysis::~TPZPostProcAnalysis(){

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
				//TPZPostProcMat * pPostProcMat = new TPZPostProcMat(matId);
				//matNumber = pcPostProcMesh->InsertMaterialObject(TPZAutoPointer<TPZMaterial>(pPostProcMat));
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
	
	pcPostProcMesh->SetDefaultOrder(pcMainMesh->GetDefaultOrder());
	pcPostProcMesh->AutoBuild();
	
	pcPostProcMesh->LoadReferred(pcMainMesh);
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

// CompEl create Functions setup

#include "pzintel.h"

#include "tpzint1point.h"
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



void TPZPostProcAnalysis::SetAllCreateFunctionsPostProc()
{
	pzgeom::TPZGeoPoint::fp = TPZPostProcAnalysis::CreatePostProcDisc;
	pzgeom::TPZGeoLinear::fp = TPZPostProcAnalysis::CreatePostProcDisc;
	pzgeom::TPZGeoQuad::fp = TPZPostProcAnalysis::CreatePostProcDisc;
	pzgeom::TPZGeoTriangle::fp = TPZPostProcAnalysis::CreatePostProcDisc;
	pzgeom::TPZGeoPrism::fp = TPZPostProcAnalysis::CreatePostProcDisc;
	pzgeom::TPZGeoTetrahedra::fp = TPZPostProcAnalysis::CreatePostProcDisc;
	pzgeom::TPZGeoPyramid::fp = TPZPostProcAnalysis::CreatePostProcDisc;
	pzgeom::TPZGeoCube::fp = TPZPostProcAnalysis::CreatePostProcDisc;
}

TPZCompEl * TPZPostProcAnalysis::CreatePostProcDisc(TPZGeoEl *gel, TPZCompMesh &mesh, int &index)
{
	return new TPZCompElPostProc< TPZCompElDisc > (mesh,gel,index);
}
