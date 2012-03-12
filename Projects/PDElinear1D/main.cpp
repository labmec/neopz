#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "pzvec.h"
#include "pzfmatrix.h"
#include "pzfstrmatrix.h"
#include "pzlog.h"
#include "pzmat1dlin.h"

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"

#include "TPZGeoLinear.h"
#include "pzgeopoint.h"

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"

#include "tpzautopointer.h"
#include "pzbndcond.h"
#include "pzanalysis.h"
#include "TPZParSkylineStructMatrix.h"
#include "pzstepsolver.h"

#include "pzstrmatrix.h"
#include "TPZFrontNonSym.h"
#include "TPZFrontSym.h"


#include <iostream>
#include <string>
#include <math.h>



//Task to DO:

//Usar TPZtransient Aanlysis
//Usar TPZtransient Material
//Do changes to contribute of Transient analysis


using namespace std;


const int matId = 1;
const int bc1 = -1;
const int bc2 = -2;

const int dirichlet = 0;
const int neumann = 1;

//
//	This program solve  a system of linear PDE - 1D I hope! (Biot poroelastic problem)
//	- a.u''(x) + b.u'(x) = f(x), em  0< x <1
//	With : u(0) = uD, du/dn = uN on x =1
//	using uD=0 e uN=0
//


// Create a Geometrical Mesh
TPZGeoMesh * MalhaGeom(int h, REAL xL, REAL  xR);

// Create a Computational Mesh
TPZCompMesh * MalhaComp(TPZGeoMesh *gmesh, int p, REAL a, REAL b);

// Assemble and Solve the generate linear System
void SolveSist(TPZAnalysis &an, TPZCompMesh *fCmesh);

// Right handside term of our Linear PDE
void ForcingFunction(const TPZVec<REAL> &pto, TPZVec<REAL> &xfloat) 
{
	double x=pto[0];
	xfloat[0] = x;
}

// Exact Solution u(x)
void SolExata(TPZVec<REAL> &pto, TPZVec<REAL> &u_exact,TPZFMatrix &du_exact)
{
	double x=pto[0];
	//	double y=pto[1];
	double xexact;
	xexact = 1.-sinh(x)/sinh(1.);
	u_exact[0] = xexact;
	//	u_exact[1] =0.;
	du_exact(0,0) = 1.-cosh(x)/sinh(1.);//dx
	//	du_exact(1,0) =0.;//dy
} 

int main(int argc, char *argv[])
{	
	std::string logs("log4cxx.doubleprojection1d");
	InitializePZLOG();
	// gRefDBase.InitializeRefPatterns();
	
	std::ofstream arg1("gmesh.txt");
	std::ofstream arg2("cmesh.txt");	
	std::ofstream file("Solution.txt",ios::app);	
	std::ofstream outerror("erros.txt",ios::app);
	
	for(int p = 1; p < 2; p++)
	{
		for(int h = 1; h <10 	; h++)
		{	
			outerror << "\n\n Solving with P = " << p << " E H = " << pow(2.,h) << "\n";
			
			TPZGeoMesh * gmesh = MalhaGeom(h, 0., 1.);
			//int p=1;
			REAL a=1.;
			REAL b=1.;
			
			TPZCompMesh * cmesh = MalhaComp(gmesh,p,a,b);
			TPZAnalysis an(cmesh);
			
			SolveSist(an,cmesh);
			
			// Print Geometrical Mesh	
			gmesh->Print(arg1);
			
			// Print Computational mesh			
			cmesh->Print(arg2);
			
			// Print Solution			
			TPZFMatrix toprint = an.Solution();
			toprint.Print("Solution", file);
			
			
			//			// Saida para mathematica
			//			std::ofstream outMath("Mathematica.txt",ios::app);
			//			std::cout << "Saida = {";
			//			for(int i = 0;  i< cmesh->ElementVec().NElements(); i++)
			//			{
			//				TPZCompEl * cel = cmesh->ElementVec()[i];
			//				TPZInterpolationSpace * sp = dynamic_cast <TPZInterpolationSpace*>(cel);
			//				TPZVec<REAL> qsi(1,0.), sol(1,0.), outfem(1,0.);
			//				TPZFMatrix axes(1,3,0.), dsol(1,1,0.);
			//				if(i != 0) std::cout << ",";
			//				for(int j = 0; j < 11; j++)
			//				{
			//					qsi[0] = -1.+2.*j/10.;
			//					sp->Solution(qsi,1,sol);
			//					gmesh->ElementVec()[i]->X(qsi,outfem);
			//					std::cout << "{" << 0 << "," << 0<<"}";
			//					if(j != 0) std::cout << ",";
			//				}
			//				if(i== cmesh->ElementVec().NElements()-1) std::cout<<"};"<<std::endl;
			//			}
			//			std::cout<<"ListPlot[Saida,Joined->True]"<<endl;
			
			
			///3. Plotar as normas dos erros 
			an.SetExact(SolExata);
			TPZVec<REAL> posproc;
			an.PostProcess(posproc,outerror); // Calcula os erros.(pzanalysis.cpp)
			outerror<<endl;
			
			//saiada para VTK
			//TPZManVector<std::string,10> scalnames(1), vecnames(1);
			//			scalnames[0] = "Solution";
			//			vecnames[0] = "Derivate";
			//			
			//			std::string plotfile("saidaT.vtk");
			//			const int dim = 2;
			//			int div = 0;
			//			an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
			//			an.PostProcess(div,dim);
			//			std::ofstream out("malha.txt");
			//			an.Print("nothing",out);
		}
	}
	
	return EXIT_SUCCESS;
}

TPZGeoMesh * MalhaGeom(int h, REAL xL, REAL  xR){
	
	TPZGeoMesh * gmesh = new TPZGeoMesh;
	int Qnodes = 2;
	gmesh->SetMaxNodeId(Qnodes-1);
	gmesh->NodeVec().Resize(Qnodes);
	TPZVec<TPZGeoNode> Node(Qnodes);
	
	TPZVec <int> TopolLine(2);
	TPZVec <int> TopolPoint(1);
	
	//indice dos nos
	int id = 0;
	int ndiv = 1;
	REAL dx = fabs(xR - xL)/ndiv;
	REAL pointco;
	for(int xi = 0; xi < Qnodes; xi++)
	{
		pointco = xL+xi*dx;
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0 ,pointco);//coord X
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}
	
	//indice dos elementos
	id=0;
	TopolPoint[0] = 0;
	new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (id,TopolPoint,bc1, *gmesh);
	id++;
	
	for (int eli=0; eli<ndiv; eli++) {
		TopolLine[0] = eli;
		TopolLine[1] = eli+1;
		new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,matId,*gmesh);
		id++;
	}
	
	TopolPoint[0] = Qnodes-1;
	new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (id,TopolPoint,bc2,*gmesh);
	
	gmesh->BuildConnectivity();
	
	
	//refinamento uniforme
	for (int ref = 0; ref < h; ref++){
		TPZVec<TPZGeoEl *> filhos;
		int n = gmesh->NElements();
		for ( int i = 0; i < n; i++ ){
			TPZGeoEl * gel = gmesh->ElementVec() [i];
			if (gel->Dimension() == 1) gel->Divide (filhos);
		}//for i
	}//ref
	
	
	//	ofstream arg("gmesh.txt");
	//	gmesh->Print(arg);
	
	return gmesh;
}

TPZCompMesh * MalhaComp(TPZGeoMesh *gmesh, int p, REAL a, REAL b){
	
	/// criar materiais
	int dim = 1;
	
	TPZMat1dLin *material;
	material = new TPZMat1dLin(matId); 
	
	TPZFMatrix xkin(1,1,1.0);
	TPZFMatrix xcin(1,1,0.0);
	TPZFMatrix xbin(1,1,1.0);
	TPZFMatrix xfin(1,1,0.0);
	
	material->SetMaterial(xkin,xcin,xbin,xfin);
	
	//	TPZMatPoisson3d *material;
	//	material = new TPZMatPoisson3d(matId,dim); 
	TPZAutoPointer<TPZMaterial> mat(material);
	
	
	
	TPZVec<REAL> convdir(3,0.);
	convdir[0]=1.;
	//REAL flux = 0.;
	
	//	material->SetParameters(diff, conv, convdir);
	material->SetForcingFunction(new TPZDummyFunction(ForcingFunction));
	
	TPZCompEl::SetgOrder(p);
	TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDimModel(dim);
	cmesh->SetAllCreateFunctionsContinuous();
	cmesh->InsertMaterialObject(mat);
	
	
	///Inserir condicao de contorno
	REAL uD=0.;
//	REAL uN=1-cosh(1.)/sinh(1.);
	TPZFMatrix val1(1,1,0.), val2(1,1,uD);
	TPZAutoPointer<TPZMaterial> BCond1 = material->CreateBC(mat, bc1,dirichlet, val1, val2);
	cmesh->InsertMaterialObject(BCond1);
	
	TPZFMatrix val12(1,1,0.), val22(1,1,uD);
	TPZAutoPointer<TPZMaterial> BCond2 = material->CreateBC(mat, bc2,dirichlet, val12, val22);
	cmesh->InsertMaterialObject(BCond2);
	
	//Ajuste da estrutura de dados computacional
	cmesh->AutoBuild();
	cmesh->AdjustBoundaryElements(); 
	cmesh->CleanUpUnconnectedNodes();
	
	//	ofstream arg("cmesh.txt");
	//	cmesh->Print(arg);
	
	return cmesh;
	
}

void SolveSist(TPZAnalysis &an, TPZCompMesh *fCmesh)
{
	
	// Caso simetrico	
	TPZSkylineStructMatrix full(fCmesh); //caso simetrico	
	an.SetStructuralMatrix(full);
	TPZStepSolver step;
	step.SetDirect(ELDLt); //caso simetrico
	an.SetSolver(step);
	an.Run();
	
	
	//	// Caso nao simetrico
	//	TPZBandStructMatrix full(fCmesh); //caso nao simetrico	
	//	an.SetStructuralMatrix(full);
	//	TPZStepSolver step;
	//	step.SetDirect(ELU); //caso nao simetrico
	//	an.SetSolver(step);
	//	an.Run();	
	
	
	//	TPZAutoPointer<TPZMaterial> mat = fCmesh->FindMaterial(matId);
	//	TPZMatPoisson3d * aximat1 = dynamic_cast<TPZMatPoisson3d*>(mat.operator->());
	//	
	//	ofstream file("Solution.out");
	//	an.Solution().Print("solution", file);    //Solution visualization on Paraview (VTK)
	//	
	
	
	//  #ifdef VTK
	//	TPZManVector<std::string,10> scalnames(1), vecnames(1);
	//	scalnames[0] = "Solution";
	//	vecnames[0] = "Derivate";
	//	
	//	std::string plotfile("saidaT.vtk");
	//	const int dim = 2;
	//	int div = 0;
	//	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	//	an.PostProcess(div,dim);
	//	std::ofstream out("malha.txt");
	//	an.Print("nothing",out);
	//  #endif	
	
	
}

