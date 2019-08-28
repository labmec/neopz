#include "pzmganalysis.h"
#include "pzcompel.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include "pztransfer.h"

#include "pzgengrid.h"
#include "pzpoisson3d.h"
#include "pzbndcond.h"

#include <fstream>

using namespace std;

TPZCompMesh *CreateMesh();

int main() {
	
	TPZCompEl::SetgOrder(1);
	TPZCompMesh *cmesh = CreateMesh();
	TPZFMatrix sol(cmesh->NEquations(),1);
	ofstream out("output.txt");
	TPZMGAnalysis mgan(cmesh);
	TPZSkylineStructMatrix strskyl(cmesh);
	mgan.SetStructuralMatrix(strskyl);
	TPZStepSolver direct;
	direct.SetDirect(ELDLt);
	mgan.SetSolver(direct);
	mgan.Run();
	int nmeshes = 2;
	TPZCompMesh *cmesh2 = 0;
	TPZGeoMesh *gmesh = 0;
	for (int imesh=0; imesh<nmeshes; imesh++) {

		if(imesh == nmeshes-1)
		{
			TPZCompEl::SetgOrder(2);
		}
		gmesh = cmesh->Reference();
		int nel = gmesh->ElementVec().NElements();
		int el;
		TPZVec<TPZGeoEl *> sub;
		for(el=0; el<nel; el++) {
			TPZGeoEl *gel = gmesh->ElementVec()[el];
			if(!gel) continue;
			gel->Divide(sub);
		}
		gmesh->ResetReference();
		cmesh2 = new TPZCompMesh(gmesh);
		cmesh->CopyMaterials(*cmesh2);
		cmesh2->AutoBuild();
		mgan.AppendMesh(cmesh2);
		mgan.Run();
		cmesh = cmesh2;

	}

	TPZTransfer trf;
	cmesh2->BuildTransferMatrix(*cmesh,trf);
	trf.Print("Transfer Matrix",out);
	TPZFMatrix sol2(cmesh2->NEquations(),1,0.);
	sol = cmesh->Solution();
	trf.TransferSolution(sol,sol2);
	cmesh2->LoadSolution(sol2);
	gmesh->Print(out);
	cmesh->Print(out);
	cmesh2->Print(out);
	cmesh->Solution().Print("Coarse mesh solution",out);
	cmesh2->Solution().Print("Fine mesh solution",out);
	
	TPZVec<REAL> ervec,truervec;
	TPZMGAnalysis::MeshError(cmesh2,cmesh,ervec,mgan.fExact,truervec);
	int i;
	cout << "TPZMGAnalysis the error between both meshes is equal to \n";
	for(i=0; i<ervec.NElements(); i++) cout << ervec[i] << ' ';
	cout << endl;
	return 0;
}

TPZCompMesh *CreateMesh()
{
	TPZGeoMesh *gmesh = new TPZGeoMesh;

	TPZManVector<int> nx(2,2);
	TPZManVector<REAL> x0(3,0.), x1(3,1.);
	x1[2]=0.;
	TPZGenGrid gen(nx,x0,x1);
	gen.Read(*gmesh);
	gen.SetBC(gmesh, 1, -1);
	TPZVec<int> corner(1,0);
	int index;
	gmesh->CreateGeoElement(EPoint, corner, -2, index);
	gmesh->BuildConnectivity();
	
	TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
	int matid = 1;
	int dimension = 2;
	TPZMatPoisson3d *poisson = new TPZMatPoisson3d(matid,dimension);
	poisson->SetInternalFlux(1.);
	TPZAutoPointer<TPZMaterial> matpoisson(poisson);
	cmesh->InsertMaterialObject(matpoisson);
	
	int bndmatid = -1;
	int bctype = 2;
	TPZFNMatrix<9> val1(1,1,1.),val2(1,1,1.);
	TPZBndCond *bndcond = new TPZBndCond(matpoisson,bndmatid,bctype,val1,val2);
	TPZAutoPointer<TPZMaterial> matbnd(bndcond);
	cmesh->InsertMaterialObject(matbnd);
	
	bndmatid = -2;
	bctype = 2;
	bndcond = new TPZBndCond(matpoisson,bndmatid,bctype,val1,val2);
	TPZAutoPointer<TPZMaterial> matbnd2(bndcond);
	cmesh->InsertMaterialObject(matbnd2);
	
	cmesh->AutoBuild();
	return cmesh;
}
