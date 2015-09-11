/**
 * @file
 * @brief A First contact with complex state-variables in . An electromagnetic wave
 * incides in a dielectric slab with a metal backing. Its
 * reflection coefficient is analised in function of the
 * angle of incidence
 * @author Francisco Orlandini
 * @since 2015
 */


#include <iostream>
#include <fstream>
#include <string>
#include "pzgmesh.h"
#include "TPZVTKGeoMesh.h"
#include "pzanalysis.h"
#include "pzbndcond.h"
#include "TPZMatValidacaoHCurlFran1.h"
#include "pzlog.h"
#include "TPZTimer.h"

//------First Validation of a HCurl Formulation-----------------


/**
 * @brief Creates gmesh corresponding to the simple domain
 * @param hDomain height of the simulation domain
 * @param wDomain width of the simulation domain
 */
TPZGeoMesh *CreateGMesh(REAL hDomain, REAL wDomain, int xDiv, int yDiv);

/**
 * @brief Creates cmesh
 * @note All the relevant parameters are arguments of this function (except constants)
 * @param gmesh the geometric mesh
 * @param pOrder polynomial approximation order
 * @param ur relative permeability of the dielectric
 * @param er relative permittivity of the dielectric
 * @param freq frequency of the plane-wave
 */
TPZCompMesh *CMesh(TPZGeoMesh *gmesh, int pOrder, STATE (& ur)( TPZVec<REAL>),STATE (& er)( TPZVec<REAL>), REAL freq);

/**
 * @brief Implements permeability of possibly inhomongeneous media
 * @param x spatial coordinates
 */
STATE urSubs(TPZVec<REAL> x);

/**
 * @brief Implements permittivity of possibly inhomongeneous media
 * @param x spatial coordinates
 */
STATE erSubs(TPZVec<REAL> x);

/**
 * @brief forcing function for testing purposes
 * @param x spatial coordinates
 * @out function value at x
 */
void forcedRHS (const TPZVec<REAL> &x, TPZVec<STATE> &out);

int main(int argc, char *argv[])
{
  TPZTimer timer;
#ifdef LOG4CXX
  InitializePZLOG();
#endif
  
  //PARAMETROS FISICOS DO PROBLEMA
  REAL freq = 1./(2 * M_PI);
  //PARAMETROS DA GEOMETRIA
  REAL hDomain = 1.;
  REAL wDomain = 1.;
  
  //PARAMETROS DE SIMULACAO
  int pOrder = 1; //ordem polinomial de aproximacao
  int dim = 2;
  int xDiv = 16;
  int yDiv = 16;
  
  timer.start();
  
	TPZGeoMesh *gmesh = CreateGMesh(hDomain, wDomain, xDiv, yDiv); //funcao para criar a malha geometrica
	
	std::ofstream out0("gmeshOriginal.vtk"); //define arquivo de saida para impressao da malha no paraview
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out0, true); //imprime a malha no formato vtk
  gmesh->Print();
  
  TPZCompMesh *cmesh = CMesh(gmesh, pOrder, urSubs , erSubs , freq); //funcao para criar a malha computacional
  bool optimizeBandwidth = true; //impede a renumeracao das equacoes do problema(para obter o mesmo resultado do Oden)
  TPZAnalysis an(cmesh,optimizeBandwidth);
  an.Run();
  
  
  an.Rhs().Print("FelPZ",std::cout,EMathematicaInput);
  timer.stop();
  
  TPZFMatrix<STATE> solucao=cmesh->Solution();//Pegando o vetor de solucao, alphaj
  solucao.Print("Sol",std::cout,EMathematicaInput);//imprime na formatacao do Mathematica
  
  //fazendo pos processamento para paraview
 TPZStack<std::string> scalnames, vecnames;
  vecnames.Push("E");//setando para imprimir campoeletrico
  std::string plotfile= "ValidacaoHCurlFran1EField.vtk";//arquivo de saida que estara na pasta debug
  an.DefineGraphMesh(dim, scalnames, vecnames, plotfile);//define malha grafica
  int postProcessResolution = 2;//define resolucao do pos processamento
  an.PostProcess(postProcessResolution);//realiza pos processamento*)


  std::cout <<"Tempo de simulacao total = "<<timer.seconds()<<" s\n";
	std::cout << "FINISHED!" << std::endl;
	
	return 0;
}



TPZGeoMesh *CreateGMesh(REAL hDomain, REAL wDomain, int xDiv, int yDiv)
{
	TPZGeoMesh * gmesh = new TPZGeoMesh;//Inicializa objeto da classe TPZGeoMesh
	
  int matId = 1; //define id para um material(formulacao fraca)
  int bc0 = -1; //define id para um material(cond contorno dirichlet)
  TPZManVector<REAL> xCoord(xDiv+1,0.0);
  TPZManVector<REAL> yCoord(yDiv+1,0.0);
  
  int nel = xDiv * yDiv;
  int nNodes = (xDiv + 1) * (yDiv + 1);
  
  //coordenadas dos nos da malha sem refinamento
  for (int i = 0; i<= xDiv; i++) {
    xCoord[i] = -1*(1-(REAL)i/xDiv)+1*((REAL)i/xDiv);
  }
  for (int i = 0; i<= yDiv; i++) {
    yCoord[i] = -1*(1-(REAL)i/yDiv)+1*((REAL)i/yDiv);
  }

  gmesh->NodeVec().Resize(nNodes); //Redimensiona o tamanho do vetor de nos da malha geometrica
  TPZVec<REAL> coord(3,0.0);
  for (long i = 0 ; i < nNodes; i++){
    coord[0] = xCoord[i%(xDiv+1)];
    coord[1] = yCoord[i/(xDiv+1)];

    gmesh->NodeVec()[i].SetCoord(coord); //seta coordenada de um no no vetor de nos da malha
    gmesh->NodeVec()[i].SetNodeId(i);
  }
#ifdef DEBUG
  for (int i = 0; i< nNodes; i++) {
    gmesh->NodeVec()[i].Print();
  }
#endif
	// Criando Elementos
  TPZVec <long> topolQuad(4); //vetor que sera inicializado com o indice dos nos de um elemento quadrado bidimensional
	TPZVec <long> topolLine(2); //vetor que sera inicializado com o indice dos nos de um elemento unidimensional
	long index; //id do elemento que sera preenchido pelo metodo CreateGeoElement
	for (long iel = 0; iel < nel; iel++) {
    const long ino1 = iel % (xDiv) + (iel / xDiv) * (xDiv + 1);
    const long ino2 = iel % (xDiv) + (iel / xDiv + 1) * (xDiv + 1);
    const long ino3 = iel % (xDiv) + (iel / xDiv + 1) * (xDiv + 1) + 1;
    const long ino4 = iel % (xDiv) + (iel / xDiv) * (xDiv + 1) + 1;
		topolQuad[0] = ino1;
		topolQuad[1] = ino2;
    topolQuad[2] = ino3;
    topolQuad[3] = ino4;
    //std::cout <<topolQuad[0]<<" "<<topolQuad[1]<<" "<<topolQuad[2]<<" "<<topolQuad[3]<<std::endl;
    gmesh->CreateGeoElement(EQuadrilateral, topolQuad, matId, index);//cria elemento quadrilateral
    gmesh->ElementVec()[index];
	}
	
  //CONDICOES DE CONTORNO
  for (int i = 0; i < 2 * xDiv; i++)
  {
    topolLine[0] = i % xDiv + (i / xDiv) * yDiv * (xDiv + 1);
    topolLine[1] = i % xDiv + (i / xDiv) * yDiv * (xDiv + 1) + 1;
    //std::cout <<topolLine[0]<<" "<<topolLine[1]<<std::endl;
    gmesh->CreateGeoElement(EOned, topolLine, bc0, index);//cria elemento unidimensional
    gmesh->ElementVec()[index];
  }
  
  for (int i = 0; i < 2 * yDiv; i++)
  {
    topolLine[0] = (i % yDiv) * (xDiv + 1) + (i / yDiv ) * xDiv;
    topolLine[1] = (i % yDiv + 1) * (xDiv + 1) + (i / yDiv ) * xDiv;
    std::cout <<topolLine[0]<<" "<<topolLine[1]<<std::endl;
    gmesh->CreateGeoElement(EOned, topolLine, bc0, index);//cria elemento unidimensional
    gmesh->ElementVec()[index];
  }
  
  
	gmesh->BuildConnectivity(); //constroi a conectividade de vizinhanca da malha
	
	return gmesh;
}


TPZCompMesh *CMesh(TPZGeoMesh *gmesh, int pOrder, STATE (& ur)( TPZVec<REAL>),STATE (& er)( TPZVec<REAL>), REAL freq)
{
  const int dim = 2; //dimensao do problema
  const int matId = 1; //define id para um material(formulacao fraca)
  const int bc0 = -1; //define id para um material(cond contorno dirichlet)
  const int dirichlet = 0, neumann = 1, mixed = 2; //tipo da condicao de contorno do problema ->default dirichlet na esquerda e na direita
  
  // Criando material
  TPZMatValidacaoHCurlFran1 *material = new TPZMatValidacaoHCurlFran1(matId,freq, ur,er);//criando material que implementa a formulacao fraca do problema de validacao
  
  material->SetForcingFunction(forcedRHS);
  ///criar malha computacional
  TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDefaultOrder(pOrder);//seta ordem polimonial de aproximacao
  cmesh->SetDimModel(dim);//seta dimensao do modelo
  // Inserindo material na malha
  cmesh->InsertMaterialObject(material);
		
  ///Inserir condicao de contorno condutores
  TPZFMatrix<STATE> val1(1,1,0.), val2(1,1,0.);
  TPZMaterial * BCond0 = material->CreateBC(material, bc0, dirichlet, val1, val2);//cria material que implementa a condicao de contorno de dirichlet
  
  
  cmesh->InsertMaterialObject(BCond0);//insere material na malha
  
  cmesh->SetAllCreateFunctionsHDiv();//define espaco de aproximacao
  
  //Cria elementos computacionais que gerenciarao o espaco de aproximacao da malha
  cmesh->AutoBuild();
  
  return cmesh;
}

STATE urSubs(TPZVec<REAL> x)
{
  STATE val(1.);
  return val;
}

STATE erSubs(TPZVec<REAL> x)
{
  STATE val(1.);
  return val;
}

void forcedRHS (const TPZVec<REAL> &x, TPZVec<STATE> &out)
{
  out[0] = -1. * ( (x[1] - 1.) * (x[1] + 1.) + 2. );
  out[1] = -1. * ( (x[0] - 1.) * (x[0] + 1.) + 2. );
  out[2] = 0;
}