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
#include "TPZMatStripline.h"
#include "TPZTimer.h"

//bibliotecas para refinamento geometrico
#include "TPZRefPatternTools.h"
#include "TPZRefPatternDataBase.h"
#include "TPZRefPattern.h"

//------Plane-wave reflection by a metal-backed dielectric slab-----------------



/**
 * @brief Creates gmesh corresponding to the stripline cavity
 * @param hDomain height of the simulation domain (conductor planes)
 * @param wDomain width of the simulation domain (absorbing boundary conditions)
 * @param hStrip height of the stripline
 * @param wStrip width of the stripline
 */
TPZGeoMesh *CreateGMesh(REAL hDomain, REAL wDomain, REAL hStrip, REAL wStrip);


/**
 * @brief Uniformly refine geomesh
 * @param gmesh geomesh to be refined
 * @param nlevel number of refinements
 */
void UniformRefineMesh(TPZGeoMesh *gmesh, int nlevel);

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

int main(int argc, char *argv[])
{
  TPZTimer timer;
  
  //PARAMETROS FISICOS DO PROBLEMA
  REAL freq = 1.0e9;
  REAL lambda=M_C/freq;
  //PARAMETROS DA GEOMETRIA
  REAL hDomain = 3.2e-3;
  REAL wDomain = 20.0e-3;
  REAL hStrip = 35.0e-6;
  REAL wStrip = 2.0e-3;
  
  //PARAMETROS DE SIMULACAO
  int pOrder = 1; //ordem polinomial de aproximacao
  int dim = 2;
  
  
  timer.start();
  
	TPZGeoMesh *gmesh = CreateGMesh(hDomain, wDomain, hStrip, wStrip); //funcao para criar a malha geometrica
	
	std::ofstream out0("gmeshOriginal.vtk"); //define arquivo de saida para impressao da malha no paraview
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out0, true); //imprime a malha no formato vtk
  gmesh->Print();
  
  UniformRefineMesh(gmesh,2);
  std::ofstream out1("gmeshUniformRefinement.vtk"); //define arquivo de saida para impressao da malha no paraview
  TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out1, true); //imprime a malha no formato vtk
  gmesh->Print();
  
  TPZCompMesh *cmesh = CMesh(gmesh, pOrder, urSubs , erSubs , freq); //funcao para criar a malha computacional
  bool optimizeBandwidth = true; //impede a renumeracao das equacoes do problema(para obter o mesmo resultado do Oden)
  TPZAnalysis an(cmesh,optimizeBandwidth);
  an.Run();
  
  timer.stop();
  
  TPZFMatrix<STATE> solucao=cmesh->Solution();//Pegando o vetor de solucao, alphaj
  solucao.Print("Sol",std::cout,EMathematicaInput);//imprime na formatacao do Mathematica
  
  //fazendo pos processamento para paraview
  TPZStack<std::string> scalnames, vecnames;
  scalnames.Push("State");//setando para imprimir u
  std::string plotfile= "StriplineEField.vtk";//arquivo de saida que estara na pasta debug
  an.DefineGraphMesh(dim, scalnames, vecnames, plotfile);//define malha grafica
  int postProcessResolution = 0;//define resolucao do pos processamento
  an.PostProcess(postProcessResolution);//realiza pos processamento


  std::cout <<"Tempo de simulacao total = "<<timer.seconds()<<" s\n";
	std::cout << "FINISHED!" << std::endl;
	
	return 0;
}



TPZGeoMesh *CreateGMesh(REAL hDomain, REAL wDomain, REAL hStrip, REAL wStrip)
{
	TPZGeoMesh * gmesh = new TPZGeoMesh;//Inicializa objeto da classe TPZGeoMesh
	
  int matId = 1; //define id para um material(formulacao fraca)
  int bc0 = -1; //define id para um material(cond contorno condutores)
  int bc1 = -2; //define id para um material(cond contorno absorcao)
  int bc2 = -3; //define id para um material(cond contorno stripline)
  
  //coordenadas dos nos da malha sem refinamento
  const REAL xCoord[8] = {-wDomain/2, wDomain/2, wDomain/2, -wDomain/2, -wStrip/2, wStrip/2, wStrip/2, -wStrip/2};
  const REAL yCoord[8] = {hDomain/2, hDomain/2, -hDomain/2, -hDomain/2, hStrip/2, hStrip/2, -hStrip/2, -hStrip/2};
  
  gmesh->NodeVec().Resize(8); //Redimensiona o tamanho do vetor de nos da malha geometrica
  TPZVec<REAL> coord(3,0.0);
  for (long i = 0 ; i < 8; i++){
    coord[0] = xCoord[i];
    coord[1] = yCoord[i];
    
    gmesh->NodeVec()[i].SetCoord(coord); //seta coordenada de um no no vetor de nos da malha
    gmesh->NodeVec()[i].SetNodeId(i);
  }
	
	// Criando Elementos
  TPZVec <long> topolQuad(4); //vetor que sera inicializado com o indice dos nos de um elemento quadrado bidimensional
	TPZVec <long> topolLine(2); //vetor que sera inicializado com o indice dos nos de um elemento unidimensional
//  TPZVec <long> TopolPoint(1); //vetor que sera inicializado com o indice do no de um elemento zero-dimensional
	long index; //id do elemento que sera preenchido pelo metodo CreateGeoElement
	for (long iel = 0; iel < 4; iel++) {
    const long ino1 = iel % 4;
    const long ino2 = (iel + 1) % 4;
    const long ino3 = (iel + 1) % 4 + 4;
    const long ino4 = iel % 4 + 4;
		topolQuad[0] = ino1;
		topolQuad[1] = ino2;
    topolQuad[2] = ino3;
    topolQuad[3] = ino4;
    //std::cout <<topolQuad[0]<<" "<<topolQuad[1]<<" "<<topolQuad[2]<<" "<<topolQuad[3]<<std::endl;
    gmesh->CreateGeoElement(EQuadrilateral, topolQuad, matId, index);//cria elemento quadrilateral
    gmesh->ElementVec()[index];
	}
	
  //CONDICOES DE CONTORNO
  int bcIds[8] = { bc0, bc1, bc0, bc1, bc2, bc2, bc2, bc2};
  for (int i = 0; i < 8; i++)
  {
    topolLine[0] = i;
    topolLine[1] = (i + 1) % 4 + 4 * (i / 4);
    //std::cout <<topolLine[0]<<" "<<topolLine[1]<<std::endl;
    gmesh->CreateGeoElement(EOned, topolLine, bcIds[i], index);//cria elemento unidimensional
    gmesh->ElementVec()[index];
  }
  
  
  
	gmesh->BuildConnectivity(); //constroi a conectividade de vizinhanca da malha
	
	return gmesh;
}

void UniformRefineMesh(TPZGeoMesh *gmesh, int nlevel)
{
  
  for (int nRefinements = 0; nRefinements < nlevel; nRefinements++) {
    long NElem = gmesh->NElements();
    TPZAdmChunkVector <TPZGeoEl * > ElemVec = gmesh->ElementVec();
    for(int i = 0; i < NElem; i++){

      TPZGeoEl *GeoEl = gmesh->Element(i);
      TPZVec<TPZGeoEl *> sons;
      GeoEl->Divide(sons);
      
    }
    
    for(int i = 0; i < gmesh->NNodes(); i++){
      
      gmesh->NodeVec()[i].SetNodeId(i);
      
    }
    
    gmesh->BuildConnectivity();
  }
  
}

TPZCompMesh *CMesh(TPZGeoMesh *gmesh, int pOrder, STATE (& ur)( TPZVec<REAL>),STATE (& er)( TPZVec<REAL>), REAL freq)
{
  const int dim = 2; //dimensao do problema
  const int matId = 1; //define id para um material(formulacao fraca)
  const int bc0 = -1; //define id para um material(cond contorno condutores)
  const int bc1 = -2; //define id para um material(cond contorno absorcao)
  const int bc2 = -3; //define id para um material(cond contorno stripline)
  const int dirichlet = 0, neumann = 1, mixed = 2; //tipo da condicao de contorno do problema ->default dirichlet na esquerda e na direita
  
  // Criando material
  TPZMatStripline *material = new TPZMatStripline(matId,freq, ur,er);//criando material que implementa a formulacao fraca do problema modelo
  
  ///criar malha computacional
  TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDefaultOrder(pOrder);//seta ordem polimonial de aproximacao
  cmesh->SetDimModel(dim);//seta dimensao do modelo
  
  // Inserindo material na malha
  cmesh->InsertMaterialObject(material);
		
  ///Inserir condicao de contorno condutores
  TPZFMatrix<STATE> val1(1,1,0.), val2(1,1,0.);
  TPZMaterial * BCond0 = material->CreateBC(material, bc0, dirichlet, val1, val2);//cria material que implementa a condicao de contorno das placas condutoras
  
  // Condicao de contorno da stripline
  //AQUIFRAN IMPLEMENTAR CONDICAO MISTA
  TPZMaterial * BCond1 = material->CreateBC(material, bc1, mixed, val1, val2);//cria material que implementa a condicao de contorno da stripline
  
  // Condicao de contorno da stripline
  TPZMaterial * BCond2 = material->CreateBC(material, bc2, dirichlet, val1, val2);//cria material que implementa a condicao de contorno da stripline
  
  cmesh->InsertMaterialObject(BCond0);//insere material na malha
  cmesh->InsertMaterialObject(BCond1);//insere material na malha
  cmesh->InsertMaterialObject(BCond2);//insere material na malha
  
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
  STATE val(4.6);
  return val;
}