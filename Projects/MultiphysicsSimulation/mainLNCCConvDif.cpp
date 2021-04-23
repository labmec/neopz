//
//  PZ
//

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzbndcond.h"
#include "TPZInterfaceEl.h"
#include "pzinterpolationspace.h"
#include "TPZLinearConvecDiff.h"

#include "tpzgeoelrefpattern.h"
#include "TPZGeoLinear.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"

#include "pzanalysis.h"
#include "pzskylstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "TPZStructMatrix.h"
#include "pzstepsolver.h"

#include "TPZGenGrid2D.h"

#include "pzlog.h"

#include <iostream>
#include <math.h>


int mainLNCC(int argc, char *argv[])
{
	///malha geometrica
  TPZGeoMesh * gmesh = new TPZGeoMesh();

  ///Criando nós
  const int nnodes = 4;
  double coord[nnodes][2] = {{0,0},{1,0},{1.,1},{0,1}};
  for(int i = 0; i < nnodes; i++) {
    int nodind = gmesh->NodeVec().AllocateNewElement();
    TPZManVector<REAL,3> nodeCoord(3);
    nodeCoord[0] = coord[i][0];
    nodeCoord[1] = coord[i][1];
    nodeCoord[2] = 0.;
    gmesh->NodeVec()[nodind].Initialize(i,nodeCoord,*gmesh);
  }

  ///Criando elementos
  const int nel = 1;
  int els[nel][4] = {{0,1,2,3}};
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
  const int nelbc = 4;
  int bcels[nelbc][3] = {{0,1,-1},{1,2,-1},{2,3,-2},{3,0,-2}};
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

	/// criar materiais
	const int dim = 2;
	const int matid = 1;
  TPZVec<REAL> ConvDir(2);
	ConvDir[0] = 1.;
	ConvDir[1] = 1.;
	const REAL difusao = 1e-12;
	const REAL F = 0.;
	const REAL SD = 1;
	TPZLinearConvecDiff * mat = new TPZLinearConvecDiff(matid,difusao, ConvDir, F, SD);

  cmesh->InsertMaterialObject(mat);
	
  ///Condições de contorno
	TPZFMatrix<STATE> val1(1,1,0.), val2(1,1,0.);
	
	val2(0,0) = 0.;
	TPZBndCond * BCondDirichletZero = mat->CreateBC(mat, -1, 0, val1, val2);//0 = Dirichlet
  cmesh->InsertMaterialObject(BCondDirichletZero);

  val2(0,0) = +1.;
  TPZBndCond * BCondDirichletUm = mat->CreateBC(mat,-2, 0, val1, val2);//0 = Dirichlet
  cmesh->InsertMaterialObject(BCondDirichletUm);
	
	const int pOrder = 2;
	cmesh->SetDefaultOrder(pOrder);
  cmesh->SetDimModel(dim);//dim = 2
	
	///cria malha H1
	cmesh->SetAllCreateFunctionsContinuous();
  	
	///Criando elementos computacionais
	cmesh->AutoBuild();
	
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
            
	///Criando structural matrix - skyline não simétrica por causa do termo convectivo
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
            
	///Exportando para Paraview
	TPZVec<std::string> scalarVars(1), vectorVars(0);
	scalarVars[0] = "Solution";
  an.DefineGraphMesh(2,scalarVars,vectorVars,"solucao.vtk");
  int resolution = 2;
	an.PostProcess(resolution);

	delete cmesh;
	delete gmesh;
	
	std::cout << "\nFIM\n";
	
  return 0;
}


