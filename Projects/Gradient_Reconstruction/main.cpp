//
//  File.cpp
//  PZ
//
//  Created by Agnaldo Farias on 7/31/12.
//  Copyright (c) 2012 LabMec-Unicamp. All rights reserved.
//

#include <iostream>

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "tpzcompmeshreferred.h"
#include "pzreferredcompel.h"
#include "pzcompel.h"
#include "pzbndcond.h"
#include "TPZInterfaceEl.h"
#include "pzinterpolationspace.h"
#include "pzpoisson3d.h"

#include "tpzgeoelrefpattern.h"
#include "TPZGeoLinear.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"

#include "pzanalysis.h"
#include "pzskylstrmatrix.h"
#include "pzstrmatrix.h"
#include "pzstepsolver.h"
#include "pzfstrmatrix.h"
#include "pzlog.h"

#include <iostream>
#include <math.h>
using namespace std;

/*
 *Projeto para validar a reconstucao do gradiente
 *Equacao: du2dx2 + du2dy2 = 0 em (0,1)x(0,1)
 *Solucao: u(x,y) = a*x + b*y
*/

REAL const coef_a = 2.;
REAL const coef_b = 4.;
TPZVec<REAL> gradU(2,0.);
TPZVec<REAL> normal_plano(2,0.);

int const matId =1;
int const bc0=-1; //em y=0
int const bc1=-2; //em x=1
int const bc2=-3; //em y=1
int const bc3=-4; //em x=0

int const dirichlet =0;
int const neumann = 1;
int const mixed = 2;

TPZGeoMesh *GMesh(int triang_elements, int nh);
TPZCompMesh *CMesh(TPZGeoMesh *gmesh, int pOrder);
void Forcingbc0(const TPZVec<REAL> &pt, TPZVec<REAL> &disp);
void Forcingbc1(const TPZVec<REAL> &pt, TPZVec<REAL> &disp);
void Forcingbc2(const TPZVec<REAL> &pt, TPZVec<REAL> &disp);
void Forcingbc3(const TPZVec<REAL> &pt, TPZVec<REAL> &disp);

void mySolve(TPZAnalysis &an, TPZCompMesh *Cmesh);
void PosProcessamento(TPZAnalysis &an, std::string plotfile,TPZFMatrix<REAL> &gradients);

void GradientReconstruction(TPZCompMesh *cmesh,TPZFMatrix<REAL> &gradients);

int main(int argc, char *argv[])
{
    gradU[0]=coef_a;
    gradU[1]=coef_b;
    normal_plano[0]=coef_a;
    normal_plano[1]=coef_b;
    
    int p = 1;
	//primeira malha
	
	// geometric mesh (initial)
	TPZGeoMesh * gmesh = GMesh(-1,1 );
    ofstream arg1("gmesh_inicial.txt");
    gmesh->Print(arg1);
    
    // First computational mesh
	TPZCompMesh * cmesh= CMesh(gmesh,p);
    ofstream arg2("cmesh_inicial.txt");
    cmesh->Print(arg2);
      
	// Solving
    TPZAnalysis an(cmesh);
    mySolve(an, cmesh);
	
	// Computing approximation of gradient
	/** 
	 * @brief Method to reconstruct a gradient after run Solve of the analysis
	 * @param cmesh Computational mesh with solution */
    TPZFMatrix<REAL> gradients;
	// Redimension of gradients storage
	gradients.Redim(cmesh->NElements(),cmesh->Dimension());
	GradientReconstruction(cmesh,gradients);
	gradients.Print();

	// Print gradient reconstructed
    string plotfile("GradientAndSolution.vtk");
	PosProcessamento(an,plotfile,gradients);

    //solucao=cmesh->Solution();
    //solucao.Print();
    
    return 0;
}

void GradientReconstruction(TPZCompMesh *cmesh,TPZFMatrix<REAL> &gradients) {
	TPZCompEl *cel;
	int i = 0, side, nneighs;
	TPZVec<REAL> normal(3,0.0);
	TPZVec<REAL> center(3,0.);
	TPZVec<REAL> solalfa(3,0.), solbeta(3,0.);
	
	REAL measure, sidemeasure;
	for(i=0;i<cmesh->NElements();i++) {
		cel = cmesh->ElementVec()[i];
		if(!cel /*|| cel->Dimension()!=cmesh->Dimension()*/) continue;
		TPZStack<TPZCompElSide> neighs;
		for(side = 0; side < cel->Reference()->NSides(); side++) {
			neighs.Resize(0);
			TPZGeoElSide gelside(cel->Reference(),side);
			if(gelside.Dimension() != cel->Dimension() - 1) continue;
			gelside.EqualorHigherCompElementList2(neighs,1,0);
			if(!neighs.NElements()) continue;
			measure = cel->VolumeOfEl();
			// for testing, we are using the solution on centroid
			TPZVec<REAL> centeralfa(3,0.);
			TPZGeoElSide gelalfa(cel->Reference(),cel->Reference()->NSides() -1);
			gelalfa.CenterPoint(centeralfa);
			cel->Solution(centeralfa,0,solalfa);
			for(nneighs=0;nneighs<neighs.NElements();nneighs++) {
				if(neighs[nneighs].Element() == cel || neighs[nneighs].Element()->Dimension() != cel->Dimension()) continue;
				sidemeasure = 1.;   //neighs[nneighs].Reference().Area();   //gelside->Area();
				neighs[nneighs].Reference().CenterPoint(center);
				neighs[nneighs].Reference().Normal(center,gelside.Element(),neighs[nneighs].Reference().Element(),normal);
               // std::cout<<"nx = " << normal[0] <<", ny = "<<normal[1]<<std::endl;
				// for testing, we are using the solution on centroid
				TPZGeoElSide gelbeta(neighs[nneighs].Reference().Element(),neighs[nneighs].Reference().Element()->NSides() -1);
				gelbeta.CenterPoint(centeralfa);
				gelbeta.Element()->Reference()->Solution(centeralfa,0,solbeta);
				
				// Incrementando no gradients para o elemento alfa o valor por vizinho
				gradients(i,0) += (sidemeasure * normal[0] *(solbeta[0] - solalfa[0]));
				gradients(i,1) += (sidemeasure * normal[1] *(solbeta[0] - solalfa[0]));
			}
		}
		gradients(i,0) /= (2*measure);
		gradients(i,1) /= (2*measure);
	}
}

void PosProcessamento(TPZAnalysis &an, std::string plotfile,TPZFMatrix<REAL> &gradients){
	TPZManVector<std::string,10> scalnames(1), vecnames(1);
	scalnames[0] = "Solution";
    vecnames[0] = "Derivate";
    
	const int dim = 2;
	int div = 0;
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an.PostProcess(div,dim);
	// Carregar gradients como soluçao na malla computacional
	// Imprimir solucao novamente - Acertar gradients dependente dos graus de liberdade no cmesh
	// cmesh->LoadSolution(gradients);
	an.PostProcess(div,dim);
	std::ofstream out("malha.txt");
	an.Print("nothing",out);
}

TPZGeoMesh *GMesh(int triang_elements, int nh){
    
    int Qnodes = 4;
	
	TPZGeoMesh * gmesh = new TPZGeoMesh;
	gmesh->SetMaxNodeId(Qnodes-1);
	gmesh->NodeVec().Resize(Qnodes);
	TPZVec<TPZGeoNode> Node(Qnodes);
	
	TPZVec <int> TopolQuad(4);
    TPZVec <int> TopolTriang(3);
	TPZVec <int> TopolLine(2);
	
	//indice dos nos
	int id = 0;
	REAL valx, dx=1.;
	for(int xi = 0; xi < Qnodes/2; xi++)
	{
		valx = xi*dx;
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0 ,valx );//coord X
		Node[id].SetCoord(1 ,0. );//coord Y
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}
	
	for(int xi = 0; xi < Qnodes/2; xi++)
	{
		valx = 1. - xi*dx;
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0 ,valx );//coord X
		Node[id].SetCoord(1 ,1. );//coord Y
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}
    
	//indice dos elementos
	id = 0;
    
    if(triang_elements==1)
    {
        TopolTriang[0] = 0;
        TopolTriang[1] = 1;
        TopolTriang[2] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (id,TopolTriang,matId,*gmesh);
        id++;
        
        TopolTriang[0] = 2;
        TopolTriang[1] = 1;
        TopolTriang[2] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (id,TopolTriang,matId,*gmesh);
        id++;
        
        TopolLine[0] = 0;
        TopolLine[1] = 1;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc0,*gmesh);
        id++;
        
        TopolLine[0] = 2;
        TopolLine[1] = 1;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc1,*gmesh);
        id++;
        
        TopolLine[0] = 3;
        TopolLine[1] = 2;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc2,*gmesh);
        id++;
        
        TopolLine[0] = 3;
        TopolLine[1] = 0;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc3,*gmesh);
    }
    else{
        TopolQuad[0] = 0;
        TopolQuad[1] = 1;
        TopolQuad[2] = 2;
        TopolQuad[3] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,matId,*gmesh);
        id++;
        
        TopolLine[0] = 0;
        TopolLine[1] = 1;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc0,*gmesh);
        id++;
        
        TopolLine[0] = 1;
        TopolLine[1] = 2;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc1,*gmesh);
        id++;
        
        TopolLine[0] = 2;
        TopolLine[1] = 3;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc2,*gmesh);
        id++;
        
        TopolLine[0] = 3;
        TopolLine[1] = 0;
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc3,*gmesh);
    }
    
	gmesh->BuildConnectivity();
    
    //#ifdef LOG4CXX
    //	if(logdata->isDebugEnabled())
    //	{
    //        std::stringstream sout;
    //        sout<<"\n\n Malha Geometrica Inicial\n ";
    //        gmesh->Print(sout);
    //        LOGPZ_DEBUG(logdata,sout.str())
    //	}
    //#endif
    
    for ( int ref = 0; ref < nh; ref++ ){
		TPZVec<TPZGeoEl *> filhos;
		int n = gmesh->NElements();
		for ( int i = 0; i < n; i++ ){
			TPZGeoEl * gel = gmesh->ElementVec() [i];
			//if (gel->Dimension() == 2) gel->Divide (filhos);
            gel->Divide (filhos);
		}//for i
	}//ref

	return gmesh;
}


TPZCompMesh *CMesh(TPZGeoMesh *gmesh, int pOrder)
{
    /// criar materiais
	int dim = 2;
	TPZMatPoisson3d *material;
	material = new TPZMatPoisson3d(matId,dim); 
	material->NStateVariables();
    
    REAL diff=1.;
    REAL conv =0.;
    TPZVec<REAL> convdir;
    convdir.Resize(dim,0.);
    convdir[0]=0.;
    material-> SetParameters(diff, conv, convdir);
    
    REAL ff=0.;
    material->SetInternalFlux(ff);
    	
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    TPZMaterial * mat(material);
    cmesh->InsertMaterialObject(mat);
    
    TPZAutoPointer<TPZFunction<STATE> > fCC0 = new TPZDummyFunction<STATE>(Forcingbc0);
    TPZAutoPointer<TPZFunction<STATE> > fCC2 = new TPZDummyFunction<STATE>(Forcingbc2);
    //TPZAutoPointer<TPZFunction<STATE> > fCC1 = new TPZDummyFunction<STATE>(Forcingbc1);
    //TPZAutoPointer<TPZFunction<STATE> > fCC3 = new TPZDummyFunction<STATE>(Forcingbc3);
    
    ///Inserir condicao de contorno
	TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
    
	TPZMaterial * BCond0 = material->CreateBC(mat, bc0,dirichlet, val1, val2);
    TPZMaterial * BCond2 = material->CreateBC(mat, bc2,dirichlet, val1, val2);

    val2(0,0)=coef_a;
    TPZMaterial * BCond1 = material->CreateBC(mat, bc1,neumann, val1, val2);
    val2(0,0)=-coef_a;
    TPZMaterial * BCond3 = material->CreateBC(mat, bc3,neumann, val1, val2);
    
    
   
    
    BCond0->SetForcingFunction(fCC0);
    BCond2->SetForcingFunction(fCC2);
   // BCond1->SetForcingFunction(fCC1);
   // BCond3->SetForcingFunction(fCC3);
    
	cmesh->SetAllCreateFunctionsContinuous();
    cmesh->InsertMaterialObject(BCond0);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    
	cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
	
	//Ajuste da estrutura de dados computacional
	cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
	cmesh->CleanUpUnconnectedNodes();
        
    //#ifdef LOG4CXX
    //	if(logdata->isDebugEnabled())
    //	{
    //        std::stringstream sout;
    //        sout<<"\n\n Malha Computacional_2 pressure\n ";
    //        cmesh->Print(sout);
    //        LOGPZ_DEBUG(logdata,sout.str());
    //	}
    //#endif
	
	return cmesh;
}


void Forcingbc0(const TPZVec<REAL> &pt, TPZVec<REAL> &disp){
	double x = pt[0];
    disp[0]= coef_a*x;
}

void Forcingbc1(const TPZVec<REAL> &pt, TPZVec<REAL> &disp){
	double y = pt[1];
    disp[0]= coef_a + coef_b*y;
}

void Forcingbc2(const TPZVec<REAL> &pt, TPZVec<REAL> &disp){
	double x = pt[0];
    disp[0]= coef_b + coef_a*x;
}

void Forcingbc3(const TPZVec<REAL> &pt, TPZVec<REAL> &disp){
    double y = pt[1];
    disp[0]= coef_b*y;
}


void mySolve(TPZAnalysis &an, TPZCompMesh *Cmesh)
{			
	//TPZBandStructMatrix full(fCmesh);
	TPZSkylineStructMatrix full(Cmesh); //caso simetrico
	an.SetStructuralMatrix(full);
	TPZStepSolver<REAL> step;
	step.SetDirect(ELDLt); //caso simetrico
	//step.SetDirect(ELU);
	an.SetSolver(step);
	an.Run();
	
	//Saida de Dados: solucao e  grafico no VT
	ofstream file("Solutout");
	an.Solution().Print("solution", file);    //Solution visualization on Paraview (VTK)
}

