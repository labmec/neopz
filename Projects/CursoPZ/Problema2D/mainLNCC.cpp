//
//  PZ
//

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzbndcond.h"
#include "TPZInterfaceEl.h"
#include "pzinterpolationspace.h"
#include "TPZMatLaplacian.h"

#include "tpzgeoelrefpattern.h"
#include "TPZGeoLinear.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"

#include "pzanalysis.h"
#include "pzskylstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzstrmatrix.h"
#include "pzstepsolver.h"

#include "pzgengrid.h"
#include "pzfunction.h"

#include "pzlog.h"

#include <iostream>
#include <math.h>

static void SolExataSteklov(const TPZVec<REAL> &loc, TPZVec<STATE> &u, TPZFMatrix<STATE> &du){
    
    const REAL n = 0;
    
    const REAL x = loc[0];
    const REAL y = loc[1];
    const REAL r = sqrt(x*x+y*y);
    const REAL t = atan2(y,x);
    const REAL sol = pow((REAL)2,0.25 + n/2.)*pow(r,0.5 + n)*cos((0.5 + n)*t);
    u[0] = sol;
    
    du(0,0) = pow((REAL)2,-0.75 + n/2.)*(1 + 2*n)*pow(pow(x,2) + pow(y,2),-0.75 + n/2.)*(x*cos((0.5 + n)*atan2(y,x)) + y*sin((0.5 + n)*atan2(y,x)));
    du(1,0) = pow((REAL)2,-0.75 + n/2.)*(1 + 2*n)*pow(pow(x,2) + pow(y,2),-0.75 + n/2.)*(y*cos((0.5 + n)*atan2(y,x)) - x*sin((0.5 + n)*atan2(y,x)));
    
}

static void NeumannEsquerda(const TPZVec<REAL> &loc, TPZVec<STATE> &result) {   ///Jorge 2017, TPZFMatrix<STATE> &du){
    REAL normal[2] = {-1,0};
    
    TPZManVector<STATE> u(1);
    TPZFMatrix<STATE> du(2,1);
    SolExataSteklov(loc,u,du);
    
    result.Resize(1);
    result[0] = du(0,0)*normal[0]+du(1,0)*normal[1];
}

static void NeumannDireita(const TPZVec<REAL> &loc, TPZVec<STATE> &result) {   ///Jorge 2017, TPZFMatrix<STATE> &du){
    REAL normal[2] = {+1,0};
    
    TPZManVector<STATE> u(1);
    TPZFMatrix<STATE> du(2,1);
    SolExataSteklov(loc,u,du);
    
    result.Resize(1);
    result[0] = du(0,0)*normal[0]+du(1,0)*normal[1];
}

static void NeumannAcima(const TPZVec<REAL> &loc, TPZVec<STATE> &result) {   ///Jorge 2017, TPZFMatrix<STATE> &du){
    REAL normal[2] = {0,+1};
    
    TPZManVector<STATE> u(1);
    TPZFMatrix<STATE> du(2,1);
    SolExataSteklov(loc,u,du);
    
    result.Resize(1);
    result[0] = du(0,0)*normal[0]+du(1,0)*normal[1];
}

void PermeabilityFunc(const TPZVec<REAL> &x, TPZVec<STATE> &val, TPZFMatrix<STATE> &perm) {
	int i, j;
	int n = perm.Rows();
	int m = perm.Cols();
	for (i = 0; i < n; i++) {
		for (j = 0; j < m; j++) {
			if (i==j || j==i-n)
				perm(i, j) = 1.;
			else
				perm(i, j) = 0.;
		}
	}
}



int main(int argc, char *argv[])
{
	///malha geometrica
  TPZGeoMesh * gmesh = new TPZGeoMesh();
  
	


  ///Criando nós
  const int nnodes = 6;
  double coord[nnodes][2] = {{-0.5,0},{0,0},{0.,0.5},{-0.5,0.5},{0.5,0},{0.5,0.5}};
  for(int i = 0; i < nnodes; i++) {
    int nodind = gmesh->NodeVec().AllocateNewElement();
    TPZManVector<REAL,3> nodeCoord(3);
    nodeCoord[0] = coord[i][0];
    nodeCoord[1] = coord[i][1];
    nodeCoord[2] = 0.;
    gmesh->NodeVec()[nodind].Initialize(i,nodeCoord,*gmesh);
  }

  ///Criando elementos
  const int nel = 2;
  int els[nel][4] = {{0,1,2,3},{1,4,5,2}};
  for(int iel = 0; iel < nel; iel++){
    TPZManVector<int64_t,4> nodind(4);
    int64_t index;
    nodind[0] = els[iel][0];
    nodind[1] = els[iel][1];
    nodind[2] = els[iel][2];
    nodind[3] = els[iel][3];
    gmesh->CreateGeoElement(EQuadrilateral,nodind,1,index);
  }

  ///Criando elementos de contorno
  const int nelbc = 6;
  int bcels[nelbc][3] = {{0,1,-3},{1,4,-2},{4,5,-4},{5,2,-6},{2,3,-6},{3,0,-5}};
  for(int iel = 0; iel < nelbc; iel++){
  	TPZManVector<int64_t,4> nodind(2);
    int64_t index;
    nodind[0] = bcels[iel][0];
    nodind[1] = bcels[iel][1];
    int matid = bcels[iel][2];
    gmesh->CreateGeoElement(EOned,nodind,matid,index);
  }

  ///Construindo conectividade da malha
	gmesh->BuildConnectivity();

	{
		std::ofstream myfile("geoMinima.txt");
		gmesh->Print(myfile);
	}	
	
			///Refinamento uniforme da malha
	
	///Inicializando padrões de refinamento uniforme
  gRefDBase.InitializeUniformRefPattern(EOned);
	gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
	const int ndiv = 4;

	for (int i = 0; i < ndiv; i++){
    int nel = gmesh->NElements();
		for (int iel = 0; iel < nel; iel++){
			TPZGeoEl * gel = gmesh->ElementVec()[iel];
			if (!gel) continue;
			if (gel->HasSubElement()) continue;//para nao dividir elementos ja divididos
      TPZVec<TPZGeoEl*> filhos;
			gel->Divide(filhos);
    }///iel
  }///i

  {
		std::ofstream myfile("geoRefinado.txt");
		gmesh->Print(myfile);
	}
  
  ///Malha computacional
  TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
	
	///Indicacao de malha DGFem. Se false, vamos criar malha H1
	const bool DGFEM = false;

	/// criar materiais
	const int dim = 2;
	const int matid = 1;
	TPZMatLaplacian * mat = new TPZMatLaplacian(matid,dim); 
	const STATE K = 1., F = 0.;
	mat->SetParameters(K,F);
	TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(PermeabilityFunc, 5);
	dummy->SetPolynomialOrder(0);
	TPZAutoPointer<TPZFunction<STATE> > func(dummy);
	mat->SetPermeabilityFunction(func);

	if(DGFEM){
		///Formulacao nao-simetria de Baumann, Oden e Babuska sem penalizacao
		mat->SetNoPenalty();
		mat->SetNonSymmetric();
	}
  cmesh->InsertMaterialObject(mat);
  
	const int pOrder = 4;
	
  ///Condições de contorno
	TPZFMatrix<STATE> val1(1,1,0.), val2(1,1,0.);
	TPZBndCond * BCondDirichletNulo = mat->CreateBC(mat, -3, 0, val1, val2);//0 = Dirichlet
  cmesh->InsertMaterialObject(BCondDirichletNulo);

  TPZBndCond * BCondNeumannZero = mat->CreateBC(mat,-2, 1, val1, val2);//1 = Neumann
  cmesh->InsertMaterialObject(BCondNeumannZero);
	
  TPZMaterial * BCondNeumannEsq = mat->CreateBC(mat, -5, 1, val1, val2);//1 = Neumann
	BCondNeumannEsq->SetForcingFunction(NeumannEsquerda,pOrder);
  cmesh->InsertMaterialObject(BCondNeumannEsq);
	
  TPZMaterial * BCondNeumannDir = mat->CreateBC(mat, -4, 1, val1, val2);//1 = Neumann
	BCondNeumannDir->SetForcingFunction(NeumannDireita,pOrder);
  cmesh->InsertMaterialObject(BCondNeumannDir);
	
  TPZMaterial * BCondNeumannAcima = mat->CreateBC(mat, -6, 1, val1, val2);//1 = Neumann
	BCondNeumannAcima->SetForcingFunction(NeumannAcima,pOrder);
  cmesh->InsertMaterialObject(BCondNeumannAcima);
    

	cmesh->SetDefaultOrder(pOrder);
  cmesh->SetDimModel(dim);//dim = 2
	
	if(DGFEM){
		///Criando malha de Galerkin descontínuo
		cmesh->SetAllCreateFunctionsDiscontinuous();
	}
	else{
    ///cria malha H1
    cmesh->SetAllCreateFunctionsContinuous();
  }
	
	///Criando elementos computacionais
	cmesh->AutoBuild();
	
	if(DGFEM){
		///Cria elementos de interface
	  TPZCreateApproximationSpace::CreateInterfaces(*cmesh);
	}
	
	{
		std::ofstream myfile("geoPosAutoBuild.txt");
		gmesh->Print(myfile);
	}
	
	{
		std::ofstream myfile("cmesh.txt");
		cmesh->Print(myfile);
	}

	///Analysis : construção do problema algébrico e inversão do sistema
	TPZAnalysis an(cmesh);
            
	///Criando structural matrix - skyline não simétrica por causa do DGFEM usado
	TPZSkylineNSymStructMatrix matriz(cmesh); 
	an.SetStructuralMatrix(matriz);
	
	///Decomposição LU
	TPZStepSolver<STATE> step;
	step.SetDirect(ELU);
	an.SetSolver(step);

	///Assemble da matriz de rigidez e vetor de carga
	an.Assemble();
	
	///Resolução do sistema
	an.Solve();
	
	///Calculando erro da aproximacao
	an.SetExact(SolExataSteklov);///definindo solucao exata do problema
	TPZVec<REAL> erro;
	std::ofstream anPostProcessFile("postprocess.txt");
	an.PostProcess(erro,anPostProcessFile);///calculando erro
	
	std::cout << "\nErro de aproximação:\n";
	std::cout << "Norma H1 = " << erro[0] << "\nNorma L2 = " << erro[1] 
	          << "\nSeminorma H1 = " << erro[2] << "\n\n";
            
	///Exportando para Paraview
	TPZVec<std::string> scalarVars(1), vectorVars(0);
	scalarVars[0] = "Solution";
  an.DefineGraphMesh(2,scalarVars,vectorVars,"solucao.vtk");
  int resolution = 0;
	an.PostProcess(resolution);

	delete cmesh;
	delete gmesh;
	
  return 0;
}


