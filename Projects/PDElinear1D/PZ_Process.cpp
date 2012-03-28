/*
 * @file
 * @brief Contains the implementation of the functions grouping the processes to: \n
 * - Generate one-dimensional geometric mesh (gmesh)
 * - Generate a computational mesh (cmesh) associated with gmesh,
 * - Assembling and solving linear system for 1D linear diff. eq. on \f$[0.0, 1.0]\f$
 */

#include "PZ_Process.h"

using namespace std;

TPZGeoMesh *GeomMesh(int h,TPZVec<int> &matId,TPZVec<int> &bc,TPZVec<REAL> &xL,TPZVec<REAL> &xR)
{	
	if(!matId.NElements() || bc.NElements() < 2)
		return NULL;
	TPZGeoMesh *gmesh = new TPZGeoMesh;
	int Qnodes = 2;
	gmesh->SetMaxNodeId(Qnodes-1);
	gmesh->NodeVec().Resize(Qnodes);
	TPZVec<TPZGeoNode> Node(Qnodes);
	
	TPZVec <int> TopolLine(2);
	TPZVec <int> TopolPoint(1);
	
	// Creating geometric nodes
	int id = 0;
	int ndiv = 1;
	REAL dx = fabs(xR[0] - xL[0])/ndiv;
	REAL pointco;
	for(int xi = 0; xi < Qnodes; xi++) {
		pointco = xL[0]+xi*dx;
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0 ,pointco);//coord X
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}

	// Creating geometric elements (point and linear)
	id=0;
	TopolPoint[0] = 0;
	new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (id,TopolPoint,bc[0], *gmesh);
	id++;
	for (int eli=0; eli<ndiv; eli++) {
		TopolLine[0] = eli;
		TopolLine[1] = eli+1;
		new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,matId[0],*gmesh);
		id++;
	}
	TopolPoint[0] = Qnodes-1;
	new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (id,TopolPoint,bc[1],*gmesh);
	
	// Building the connectivity between elements as (element, side)
	gmesh->BuildConnectivity();
	
	// If necessary it refines level by level initial mesh
	UniformRefine(h,gmesh);
	
	return gmesh;
}

void UniformRefine(int h,TPZGeoMesh *gmesh) {
	//uniform refinement
	TPZVec<TPZGeoEl *> filhos;
	for (int ref = 0; ref < h; ref++) {
		int n = gmesh->NElements();
		for(int i = 0; i < n; i++ ) {
			TPZGeoEl *gel = gmesh->ElementVec()[i];
			gel->Divide(filhos);
		}
	}
}

TPZCompMesh *CompMesh(TPZGeoMesh *gmesh,int p, TPZMaterial *material,TPZVec<int> &bc,TPZVec<int> &bcType) {

	if(!material || bc.NElements()<2 || bcType.NElements() != bc.NElements()) return NULL;
	int dim = 1;
	
	
	TPZAutoPointer<TPZMaterial> mat(material);
	
	// related to interpolation space
	TPZCompEl::SetgOrder(p);
	TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDimModel(dim);
	cmesh->SetAllCreateFunctionsContinuous();
	cmesh->InsertMaterialObject(mat);
	
	// Related to boundary conditions
	//	REAL uN=1-cosh(1.)/sinh(1.);
	TPZFMatrix<REAL> val1(1,1,0.), val2(1,1,0.);
	if(!bcType[0])  // dirichlet
		val2.PutVal(0,0,0.0);
	TPZAutoPointer<TPZMaterial> BCond1 = material->CreateBC(mat, bc[0],bcType[0], val1, val2);
	cmesh->InsertMaterialObject(BCond1);
	
	if(!bcType[1])  // dirichlet
		val2.PutVal(0,0,0.0);
	TPZAutoPointer<TPZMaterial> BCond2 = material->CreateBC(mat, bc[1],bcType[1], val1, val2);
	cmesh->InsertMaterialObject(BCond2);

	//Adjusting data
	cmesh->AutoBuild();
	cmesh->AdjustBoundaryElements(); 
	cmesh->CleanUpUnconnectedNodes();
	
	return cmesh;
}

void SolveSist(TPZAnalysis &an, TPZCompMesh *fCmesh)
{
	// Symmetric case	
	TPZSkylineStructMatrix full(fCmesh);
	an.SetStructuralMatrix(full);
	TPZStepSolver<REAL> step;
	step.SetDirect(ELDLt);
	an.SetSolver(step);
	an.Run();

	//	Nonsymmetric case
	//	TPZBandStructMatrix full(fCmesh);
	//	an.SetStructuralMatrix(full);
	//	TPZStepSolver step;
	//	step.SetDirect(ELU);
	//	an.SetSolver(step);
	//	an.Run();

	//	TPZAutoPointer<TPZMaterial> mat = fCmesh->FindMaterial(matId);
	//	TPZMatPoisson3d * aximat1 = dynamic_cast<TPZMatPoisson3d*>(mat.operator->());
	//	
	//	ofstream file("Solution.out");
	//	an.Solution().Print("solution", file);    //Solution visualization on Paraview (VTK)
	//	
}

// Output as Mathematica format
void OutputMathematica(std::ofstream &outMath, TPZCompMesh *cmesh) {
	outMath << "Saida = {";
	int nelem = cmesh->ElementVec().NElements();
	int pointsInElement = 11;
	for(int i = 0;  i< nelem; i++)
	{
		TPZCompEl * cel = cmesh->ElementVec()[i];
		if(cel->Reference()->Dimension() < cmesh->Dimension()) continue;
		TPZInterpolationSpace * sp = dynamic_cast <TPZInterpolationSpace*>(cel);
		TPZVec<REAL> qsi(1,0.), sol(1,0.), outfem(3,0.);
		TPZFMatrix<REAL> axes(1,3,0.), dsol(1,1,0.);
		for(int j = 0; j < pointsInElement; j++)
		{
			qsi[0] = -1.+2.*j/10.;
			sp->Solution(qsi,0,sol);
			cel->Reference()->X(qsi,outfem);
			outMath << "{" << outfem[0] << "," << sol[0] << "}";
			if(j != pointsInElement-1) outMath << ",";
		}
		if(i == nelem-1) outMath << "};" << std::endl;
		else outMath << ",";
	}
	outMath << "ListPlot[Saida,Joined->True]"<< endl;
}

// Output as VTK (Visualization Tool Kit) format
void OutputVTK(std::string &outVTK, TPZCompMesh *cmesh,TPZAnalysis &an) {
	TPZManVector<std::string,10> scalnames(1), vecnames(1);
	scalnames[0] = "Solution";
	vecnames[0] = "Derivate";
	
	const int dim = 2;
	int div = 0;
	an.DefineGraphMesh(dim,scalnames,vecnames,outVTK);
	an.PostProcess(div,dim);
	std::ofstream out("malha.txt");
	an.Print("nothing",out);
}
