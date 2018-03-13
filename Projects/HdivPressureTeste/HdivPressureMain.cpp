/*
 * @file
 * @author Denise de Siqueira
 * @since 08/05/12.
 */

#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include <iostream>
#include <cstdlib>
#include "pzgengrid.h"
#include "pzgmesh.h"
#include "pzgeoelbc.h"
#include "pzcmesh.h"
#include "tpzcompmeshreferred.h"
#include "pzcompel.h"
#include "pzpoisson3d.h"
#include "pzbndcond.h"
#include "pzanalysiserror.h"
#include "pzanalysis.h"
#include "pzcmesh.h"
#include "pzstepsolver.h"
#include "TPZParFrontStructMatrix.h"
#include "pzmatrix.h"
#include "TPZCompElDisc.h"
#include "pzfstrmatrix.h"
#include "pzinterpolationspace.h"
#include "pzsubcmesh.h"
#include "pzlog.h"
#include "pzelctemp.h"
#include "pzelchdiv.h"
#include "pzshapequad.h"
#include "pzshapetriang.h"
#include "pzgeoquad.h"
#include "pzgeotriangle.h"
#include "pzfstrmatrix.h"
#include "pzgengrid.h"
#include "pzbndcond.h"
#include "TPZMaterial.h"
#include "tpzquadrilateral.h"
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "pzlog.h"
#include <cmath>

#include "TPZRefPattern.h"
#include "pzhdivpressure.h"

#ifdef LOG4CXX

static LoggerPtr logger(Logger::getLogger("pz.HdivPressure.main"));

#endif

void SolveLU ( TPZAnalysis &an );
TPZGeoMesh * MalhaGeoT(const int h);
TPZGeoMesh * MalhaGeoQ(const int h);
TPZCompMesh *CreateCompMesh2d(TPZGeoMesh &gmesh,int porder);


void Forcing1(const TPZVec<REAL> &pt, TPZVec<REAL> &disp) {
	//double x = pt[0];
//	double y = pt[1];
		disp[0]= 0.;
		return;
}
void SolExata(const TPZVec<REAL> &pt, TPZVec<REAL> &p, TPZFMatrix<REAL> &flux ) {
	double x = pt[0];
	//double y = pt[1];
	TPZVec<REAL> disp;
    p[0]= 0.5*(1-x*x);//Solucao
	flux(0,0)= (-1)*x;
	flux(1,0)=  0.;
  flux(2,0)=0.;//divergente
	
	
}
void CC1(const TPZVec<REAL> &pt, TPZVec<REAL> &f) {
		f[0] = 0.;
	
}

/** Resolver o problema do tipo 
 * -Laplac(u) = 0
 * u = u0 em todo contorno
 */

using namespace std;
int main()
{
#ifdef LOG4CXX
		InitializePZLOG();
#endif
		
		
#ifdef LOG4CXX
   // if (logger->isDebugEnabled())
		//{
				std::stringstream sout;
				sout<< "Teste de HdivPressure"<< std::endl;
				LOGPZ_DEBUG(logger, sout.str().c_str());
				// LOGPZ_DEBUG(logger,"Teste de HdivPressure")
	//	}
#endif
		for (int porder=2; porder<3; porder++) {
		
		for(int Nref=0;Nref<1;Nref++){
			
			TPZGeoMesh *geomesh = MalhaGeoQ(Nref);
			
			TPZCompMesh *cmesh = CreateCompMesh2d(*geomesh,porder);
			
			cmesh->LoadReferences();//mapeia para a malha geometrica lo
			
			TPZAdmChunkVector<TPZCompEl *> elvec = cmesh->ElementVec();
				


            
			//2. Resolve o problema
				
			TPZAnalysis analysis(cmesh);
			cmesh->SaddlePermute();
			SolveLU ( analysis );
			
						
			/*
			
			
			//4. visualizacao grafica usando vtk
			 TPZVec<std::string> scalnames(1), vecnames(1);
			 
			// scalnames[0] = "Divergence";
			// scalnames[0] = "Pressure";
			 scalnames[0] = "ExactPressure";
			 //	scalnames[2] = "ExactDiv";
			 
			 
			 vecnames[0] = "ExactFlux";
			// vecnames[1] = "Flux";
			 //scalnames[2] = "Divergence";
			 
			 
			 //vecnames[0] = "Derivate";
			 
			 
			 std::string plotfile("SolGraph.vtk");
			 const int dim = 2;
			 int div = 3;
			 analysis.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
			 analysis.PostProcess(div);
			*/
			
		}}
	
	
	return 0;
}



TPZCompMesh *CreateCompMesh2d(TPZGeoMesh &gmesh,int porder){
	TPZCompEl::SetgOrder(porder);
	TPZCompMesh *comp = new TPZCompMesh(&gmesh);
	
	// Criar e inserir os materiais na malha
	TPZMatPoisson3d *mat = new TPZMatPoisson3d(1,2);
	TPZMaterial * automat(mat);
	comp->InsertMaterialObject(automat);
	
  TPZAutoPointer<TPZFunction<STATE> > force1 = new TPZDummyFunction<STATE>(Forcing1);
	mat->SetForcingFunction(force1);
	TPZAutoPointer<TPZFunction<STATE> > exata1 = new TPZDummyFunction<STATE>(SolExata);
	mat->SetForcingFunctionExact(exata1);
	///Inserir condicoes de contorno
	
	TPZFMatrix<REAL> val1(1,1,0.),val2(1,1,0.);
	TPZFMatrix<REAL> val11(1,1,0.), val22(1,1,0.);
	TPZMaterial *bnd = automat->CreateBC (automat,-1,0,val1,val2);//1
	TPZMaterial *bnd2 = automat->CreateBC (automat,-2,0,val1,val2);
	TPZMaterial *bnd3 = automat->CreateBC (automat,-3,0,val1,val2);//1
	TPZMaterial *bnd4 = automat->CreateBC (automat,-4,0,val1,val2);
	
	TPZAutoPointer<TPZFunction<STATE> > fCC1 = new TPZDummyFunction<STATE>(CC1);
	bnd->SetForcingFunction(fCC1);
	bnd2->SetForcingFunction(fCC1);
	bnd3->SetForcingFunction(fCC1);
	bnd4->SetForcingFunction(fCC1);
	
	comp->InsertMaterialObject(bnd);
	comp->InsertMaterialObject(bnd2);
	comp->InsertMaterialObject(bnd3);
	comp->InsertMaterialObject(bnd4);	
	//espaco de aproximacao
	//comp->SetAllCreateFunctionsHDiv();
		comp->SetAllCreateFunctionsHDivPressure();
	//comp->SetAllCreateFunctionsContinuous();
	
	// Ajuste da estrutura de dados computacional
	comp->AutoBuild();
	comp->AdjustBoundaryElements();//ajusta as condicoes de contorno
	comp->CleanUpUnconnectedNodes();//deleta os nos que nao tem elemntos conectados
	comp->SetName("Malha Computacional Apos Criada");
	
#ifdef LOG4CXX
	{
		std::stringstream sout;
		comp->Print(sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	
	
	
    return comp;
	
}


TPZGeoMesh * MalhaGeoT(const int h){//malha triangulo
	
	TPZGeoMesh *gmesh = new TPZGeoMesh();
	
	//Criar ns
	const int nnode = 4;//AQUI
	const int nelem = 2;
	TPZGeoEl *elvec[nelem];	
	const int dim = 2;//AQUI
	
	REAL co[nnode][dim] = {{-1.,-1.},{1.,-1.},{1.,1.},{-1.,1.}};
	int indices[2][nnode];//como serao enumerados os nos
	
	//el 1
	indices[0][0] = 0;
	indices[0][1] = 1;
	indices[0][2] = 3;
	//el2
	indices[1][0] = 2;
	indices[1][1] = 3;
	indices[1][2] = 1;
	
	
	int nod;
	TPZVec<REAL> coord(dim);
	for(nod=0; nod<nnode; nod++) {
		int nodind = gmesh->NodeVec().AllocateNewElement();
		
		for(int d = 0; d < dim; d++)
		{
			coord[d] = co[nod][d];
		}
		gmesh->NodeVec()[nod].Initialize(nodind,coord,*gmesh);
	}
	//Criacao de elementos
	
	
	TPZVec<int> nodind1(3);
	TPZVec<int> nodind2(3);
	for(int i=0; i<3; i++){
		nodind1[i] = indices[0][i];
		nodind2[i] = indices[1][i];
	}
	
	int index;
	elvec[0] = gmesh->CreateGeoElement(ETriangle,nodind1,1,index); //AQUI
	elvec[1] = gmesh->CreateGeoElement(ETriangle,nodind2,1,index); //AQUI
	
	
	gmesh->BuildConnectivity();
	
	
	//Cria as condicoes de contorno
	TPZGeoElBC gbc1(elvec[0],3,-1);// condicao de fronteira tipo -1: 
	TPZGeoElBC gbc2(elvec[0],5,-2);// condicao de fronteira tipo -2: 
	
	TPZGeoElBC gbc3(elvec[1],3,-3);// condicao de fronteira tipo -3: 
	TPZGeoElBC gbc4(elvec[1],5,-4);// condicao de fronteira tipo -4: 
	
	const std::string nameref;
	
	TPZAutoPointer<TPZRefPattern> ref;
	
	
	//	Refinamento uniforme
	for(int ref = 0; ref < h; ref++){// h indica o numero de refinamentos
		TPZVec<TPZGeoEl *> filhos;
		int n = gmesh->NElements();
		for(int i = 0; i < n; i++){
			TPZGeoEl * gel = gmesh->ElementVec()[i];
			if(!gel->HasSubElement())
			{
				gel->Divide(filhos);
			}		
		}}
	
	
#ifdef LOG4CXX
	{
		std::stringstream sout;
		gmesh->Print(sout);
		LOGPZ_DEBUG(logger, sout.str().c_str());
	}
#endif 		 		 
	return gmesh;
}
TPZGeoMesh * MalhaGeoQ( const int h )/*QUADRILATEROS*/ 
{
	TPZGeoMesh *gmesh = new TPZGeoMesh();
		REAL co[4][2] = {{-1.,-1.},{1.,-1.},{1.,1.},{-1.,1.}};//{{0.,0.},{1.,0.},{1.,1.},{0.,1.}};
	int indices[1][4] = {{0,1,2,3}};
	
	int nnode = 4;
	const int nelem = 1;
	TPZGeoEl *elvec[nelem];
	int nod;
	for ( nod=0; nod<nnode; nod++ )
	{
		int nodind = gmesh->NodeVec().AllocateNewElement();
		TPZVec<REAL> coord ( 2 );
		coord[0] = co[nod][0];
		coord[1] = co[nod][1];
		gmesh->NodeVec() [nodind].Initialize ( nod,coord,*gmesh );
	}
	
	int el;
	for ( el=0; el<nelem; el++ )
	{
		TPZVec<int> nodind ( 4 );
		for ( nod=0; nod<4; nod++ ) nodind[nod]=indices[el][nod];
		int index;
		elvec[el] = gmesh->CreateGeoElement ( EQuadrilateral,nodind,1,index );
	}
	
	gmesh->BuildConnectivity();
	
	TPZGeoElBC gbc1 ( elvec[0],4,-1 );// condicao de fronteira tipo -1: (x,y=0)
	TPZGeoElBC gbc2 ( elvec[0],5,-2 );// condicao de fronteira tipo -2: (x=1,y)
	TPZGeoElBC gbc3 ( elvec[0],6,-3 );// condicao de fronteira tipo -3: (x,y=1)
	TPZGeoElBC gbc4 ( elvec[0],7,-4 );// condicao de fronteira tipo -4: (x=0,y)
	
	///Refinamento uniforme
	for ( int ref = 0; ref < h; ref++ )
	{// h indica o numero de refinamentos
		TPZVec<TPZGeoEl *> filhos;
		int n = gmesh->NElements();
		for ( int i = 0; i < n; i++ )
		{
			TPZGeoEl * gel = gmesh->ElementVec() [i];
			//if ( gel->Dimension() == 2 ) gel->Divide ( filhos );
			if(!gel->HasSubElement())
			{
				gel->Divide(filhos);
			}	
		}//for i
	}//ref
	
	return gmesh;
}



void SolveLU ( TPZAnalysis &an ){
	// Com matriz mal condicionada a solução é poluída com resíduos,
	// tanto com LU quanto Choleski. Isso resulta em não simetrias.
	TPZCompMesh *malha = an.Mesh();
	//TPZFrontStructMatrix<TPZFrontNonSym> mat ( malha );// não funciona com método iterativo
	TPZFStructMatrix mat( malha );
	//	TPZSpStructMatrix mat( malha );
	TPZStepSolver<REAL> solv;
	
	solv.SetDirect ( ELU );
	//		solv.SetDirect(ECholesky);
	
	cout << "ELU " << endl;
	an.SetSolver ( solv );
	an.SetStructuralMatrix ( mat );
	cout << endl;
	an.Solution().Redim ( 0,0 );
	cout << "Assemble " << endl;
	an.Assemble();
	an.Solve();
	cout << endl;
	cout << "No equacoes = " << malha->NEquations() << endl;
	// cout << "Banda = " << malha->BandWidth() << endl;
}
