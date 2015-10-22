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
#include "pzgengrid.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include "tpzhierarquicalgrid.h"
#include "pzgeoquad.h"
#include "pzgeopoint.h"
#include "tpzgeoelrefpattern.h"
//------First Validation of a HCurl Formulation-----------------
TPZFNMatrix< 9 , REAL > rotationMatrix(3,3,0.);

void CreateGMesh(TPZGeoMesh * &gmsh, int type, REAL hDomain, REAL wDomain, REAL x0, REAL y0, int xDiv, int yDiv);

/**
 * @brief Creates gmesh corresponding to the simple domain
 * @param hDomain height of the simulation domain
 * @param wDomain width of the simulation domain
 */
TPZGeoMesh *CreateRectangularGMesh(REAL hDomain, REAL wDomain, REAL x0, REAL y0, int xDiv, int yDiv);

TPZGeoMesh *CreateZigZagGMesh(REAL hDomain, REAL wDomain, REAL x0, REAL y0, int xDiv, int zDiv);



TPZGeoMesh *CreateTriangularGMesh(REAL hDomain, REAL wDomain, REAL x0, REAL y0,  int xDiv, int yDiv);

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
void createRotationMatrix (TPZFNMatrix<9,REAL> &mat, const REAL thetaRotation, const int axis);

int main(int argc, char *argv[])
{
  TPZTimer timer;
#ifdef LOG4CXX
  InitializePZLOG();
#endif
  
  rotationMatrix.Identity();
  const double thetaRot = M_PI_4;
  createRotationMatrix(rotationMatrix, thetaRot, 1);
  
  //PARAMETROS FISICOS DO PROBLEMA
  REAL freq = 1. / (2*M_PI*sqrt(M_UZERO*M_EZERO) );
  //PARAMETROS DA GEOMETRIA
  REAL hDomain = 2.;
  REAL wDomain = 2.;
  const REAL x0 = 0., y0 = 0.;
  
  //PARAMETROS DE SIMULACAO
  int pOrder = 1; //ordem polinomial de aproximacao
  int dim = 2;
  int xDiv = 16;
  int yDiv = 16;
  
  rotationMatrix.Identity();
  const bool rotating = true;
  if ( rotating)
  {
    const double thetaRot = M_PI_4;
    enum rotAxis {xAxis=0, yAxis, zAxis};
    createRotationMatrix(rotationMatrix, thetaRot, zAxis);
  }
  
  timer.start();
  enum meshType{ createRectangular, createTriangular, createZigZag};
  const int mesh = createTriangular;
  
  TPZGeoMesh *gmesh = new TPZGeoMesh();
  CreateGMesh(gmesh, mesh, hDomain, wDomain, x0, y0, xDiv, yDiv); 
  
  TPZCompMesh *cmesh = CMesh(gmesh, pOrder, urSubs , erSubs , freq); //funcao para criar a malha computacional
  bool optimizeBandwidth = true; 
  TPZAnalysis an(cmesh,optimizeBandwidth);
  //configuracoes do objeto de analise
  TPZSkylineStructMatrix skylstr(cmesh); //caso simetrico
  skylstr.SetNumThreads(0);
  //    TPZSkylineNSymStructMatrix full(fCmesh);
  an.SetStructuralMatrix(skylstr);
  TPZStepSolver<STATE> step;
  step.SetDirect(ELDLt); //caso simetrico
  an.SetSolver(step);
  TPZStack<std::string> scalnames, vecnames;
  vecnames.Push("realE");//setando para imprimir campoeletrico
  std::string plotfile= "../ValidacaoHCurlFran1EField.vtk";//arquivo de saida que estara na pasta debug
  an.DefineGraphMesh(dim, scalnames, vecnames, plotfile);//define malha grafica
  int postProcessResolution = 2 ;//define resolucao do pos processamento
  //fim das configuracoes do objeto de analise
  
  an.Run();
  
  timer.stop();
#ifdef PZDEBUG
  TPZFMatrix<STATE> solucao=cmesh->Solution();//Pegando o vetor de solucao, alphaj
  solucao.Print("Sol",std::cout,EMathematicaInput);//imprime na formatacao do Mathematica
#endif
  //AQUIFRAN APAGAR
//  TPZFNMatrix<200,STATE> sds(solucao);
//  for (int i = 0; i < sds.Rows(); i++) {
//    sds.Zero();
//    sds(i,0) = 1.;
//    cmesh->Solution() = sds;
//    an.PostProcess(postProcessResolution);//realiza pos processamento*)
//  }
  an.PostProcess(postProcessResolution);//realiza pos processamento*)


  std::cout <<"Tempo de simulacao total = "<<timer.seconds()<<" s\n";
	std::cout << "FINISHED!" << std::endl;
	
	return 0;
}

void CreateGMesh(TPZGeoMesh * &gmesh, int type, REAL hDomain, REAL wDomain, REAL x0, REAL y0, int xDiv, int yDiv)
{
  enum meshType{ createRectangular, createTriangular, createZigZag};
  
  switch (type) {
    case createRectangular:
    {
      gmesh = CreateRectangularGMesh(hDomain, wDomain, x0, y0, xDiv, yDiv); //funcao para criar a malha geometrica
      for (int i = 0; i < gmesh->NNodes() ; i ++)
      {
        TPZGeoNode node = gmesh->NodeVec()[i];
        TPZManVector<REAL,3> xRotated(3,0.), x(3,0.);
        node.GetCoordinates(x);
        xRotated[0] = rotationMatrix(0,0) * x[0] + rotationMatrix(0,1) * x[1] + rotationMatrix(0,2) * x[2];
        xRotated[1] = rotationMatrix(1,0) * x[0] + rotationMatrix(1,1) * x[1] + rotationMatrix(1,2) * x[2];
        xRotated[2] = rotationMatrix(2,0) * x[0] + rotationMatrix(2,1) * x[1] + rotationMatrix(2,2) * x[2];
        node.SetCoord(xRotated);
        gmesh->NodeVec()[i] = node;
      }
      std::ofstream outTxt("../gmeshRegular.txt"); //define arquivo de saida para impressao da malha no
      gmesh->Print(outTxt);
      std::ofstream out0("../gmeshOriginal.vtk"); //define arquivo de saida para impressao da malha no paraview
      TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out0, true); //imprime a malha no formato vtk
    }
      break;
    case createTriangular:
    {
      gmesh = CreateTriangularGMesh(hDomain, wDomain, x0, y0, xDiv, yDiv); //funcao para criar a malha geometrica
      for (int i = 0; i < gmesh->NNodes() ; i ++)
      {
        TPZGeoNode node = gmesh->NodeVec()[i];
        TPZManVector<REAL,3> xRotated(3,0.), x(3,0.);
        node.GetCoordinates(x);
        xRotated[0] = rotationMatrix(0,0) * x[0] + rotationMatrix(0,1) * x[1] + rotationMatrix(0,2) * x[2];
        xRotated[1] = rotationMatrix(1,0) * x[0] + rotationMatrix(1,1) * x[1] + rotationMatrix(1,2) * x[2];
        xRotated[2] = rotationMatrix(2,0) * x[0] + rotationMatrix(2,1) * x[1] + rotationMatrix(2,2) * x[2];
        node.SetCoord(xRotated);
        gmesh->NodeVec()[i] = node;
      }
      std::ofstream outTxt("../gmeshTriangular.txt"); //define arquivo de saida para impressao da malha no
      gmesh->Print(outTxt);
      std::ofstream out0("../gmeshTriangular.vtk"); //define arquivo de saida para impressao da malha no paraview
      TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out0, true); //imprime a malha no formato vtk
      
    }
      break;
    case createZigZag:
    {
      gmesh = CreateZigZagGMesh(hDomain, wDomain, x0, y0, xDiv, yDiv);
      for (int i = 0; i < gmesh->NNodes() ; i ++)
      {
        TPZGeoNode node = gmesh->NodeVec()[i];
        TPZManVector<REAL,3> xRotated(3,0.), x(3,0.);
        node.GetCoordinates(x);
        xRotated[0] = rotationMatrix(0,0) * x[0] + rotationMatrix(0,1) * x[1] + rotationMatrix(0,2) * x[2];
        xRotated[1] = rotationMatrix(1,0) * x[0] + rotationMatrix(1,1) * x[1] + rotationMatrix(1,2) * x[2];
        xRotated[2] = rotationMatrix(2,0) * x[0] + rotationMatrix(2,1) * x[1] + rotationMatrix(2,2) * x[2];
        node.SetCoord(xRotated);
        gmesh->NodeVec()[i] = node;
      }
      std::ofstream outTxt("../gmeshZigZag.txt"); //define arquivo de saida para impressao da malha no
      gmesh->Print(outTxt);
      std::ofstream out0("../gmeshZigZag.vtk"); //define arquivo de saida para impressao da malha no paraview
      TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out0, true); //imprime a malha no formato vtk
    }
      break;
    default:
      DebugStop();
      break;
  }
  
}

TPZGeoMesh *CreateTriangularGMesh(REAL hDomain, REAL wDomain, REAL x0, REAL y0, int xDiv, int yDiv)
{
  
  TPZManVector<int,3> nx(3,0);
  TPZManVector<REAL,3> llCoord(3,0.) , ulCoord(3,0.) , urCoord(3,0.) , lrCoord(3,0.);
  llCoord[0] = x0 - wDomain/2;
  llCoord[1] = y0 - hDomain/2;
  
  ulCoord[0] = x0 - wDomain/2;
  ulCoord[1] = y0 + hDomain/2;
  
  urCoord[0] = x0 + wDomain/2;
  urCoord[1] = y0 + hDomain/2;
  
  lrCoord[0] = x0 + wDomain/2;
  lrCoord[1] = y0 - hDomain/2;
  
  nx[0]=xDiv;
  nx[1]=yDiv;
  int numl = 1;
  REAL rot = 0.5;
  TPZGenGrid gengrid = TPZGenGrid(nx, llCoord, urCoord,  numl, rot);
  gengrid.SetElementType(ETriangle);
  
  TPZGeoMesh *gmesh = new TPZGeoMesh();
  const int matId = 1; //define id para um material(formulacao fraca)
  const int bc0 = -1; //define id para um material(cond contorno dirichlet)
  const int bc1 = -2; //define id para um material(cond contorno mista)
  
  gengrid.Read(gmesh , matId);
  
  gengrid.SetBC(gmesh, ulCoord, llCoord, bc0);
  gengrid.SetBC(gmesh, urCoord, ulCoord, bc0);
  gengrid.SetBC(gmesh, lrCoord, urCoord, bc1);
  gengrid.SetBC(gmesh, llCoord, lrCoord, bc0);
  
  return gmesh;
}

TPZGeoMesh *CreateZigZagGMesh(REAL hDomain, REAL wDomain, REAL x0, REAL y0, int xDiv, int zDiv)
{
  TPZManVector<int,3> nx(3,0);
  TPZManVector<REAL,3> llCoord(3,0.) , ulCoord(3,0.) , urCoord(3,0.) , lrCoord(3,0.);
  llCoord[0] = x0 - wDomain/2;
  llCoord[1] = y0 - hDomain/2;
  
  ulCoord[0] = x0 - wDomain/2;
  ulCoord[1] = y0 + hDomain/2;
  
  urCoord[0] = x0 + wDomain/2;
  urCoord[1] = y0 + hDomain/2;
  
  lrCoord[0] = x0 + wDomain/2;
  lrCoord[1] = y0 - hDomain/2;
  
  nx[0]=xDiv;
  nx[1]=zDiv;
  int numl = 1;
  REAL rot = 0.5;
  TPZGenGrid gengrid = TPZGenGrid(nx, llCoord, urCoord,  numl, rot);
  gengrid.SetZigZagPattern();
  gengrid.SetElementType(EQuadrilateral);
  
  TPZGeoMesh *gmesh = new TPZGeoMesh();
  const int matId = 1; //define id para um material(formulacao fraca)
  const int bc0 = -1; //define id para um material(cond contorno dirichlet)
  const int bc1 = -2; //define id para um material(cond contorno mista)
  
  gengrid.Read(gmesh , matId);
  
  gengrid.SetBC(gmesh, ulCoord, llCoord, bc0);
  gengrid.SetBC(gmesh, urCoord, ulCoord, bc0);
  gengrid.SetBC(gmesh, lrCoord, urCoord, bc1);
  gengrid.SetBC(gmesh, llCoord, lrCoord, bc0);
  
  
  return gmesh;
}

TPZGeoMesh *CreateRectangularGMesh(REAL hDomain, REAL wDomain, REAL x0, REAL y0, int xDiv, int yDiv)
{
	TPZGeoMesh * gmesh = new TPZGeoMesh;//Inicializa objeto da classe TPZGeoMesh
	
  const int matId = 1; //define id para um material(formulacao fraca)
  const int bc0 = -1; //define id para um material(cond contorno baixo)
  const int bc1 = -2; //define id para um material(cond contorno cima)
  const int bc2 = -3; //define id para um material(cond contorno esquerda)
  const int bc3 = -4; //define id para um material(cond contorno direita)
  TPZManVector<REAL> xCoord(xDiv+1,0.0);
  TPZManVector<REAL> yCoord(yDiv+1,0.0);
  
  int nel = xDiv * yDiv;
  int nNodes = (xDiv + 1) * (yDiv + 1);
  
  
  for (int i = 0; i<= xDiv; i++) {
    xCoord[i] = (x0 - wDomain/2)*(1-(REAL)i/xDiv) + (x0 + wDomain/2)*((REAL)i/xDiv);
  }
  for (int i = 0; i<= yDiv; i++) {
    yCoord[i] = (y0 - hDomain/2)*(1-(REAL)i/yDiv) + (y0 + hDomain/2)*((REAL)i/yDiv);
  }

  gmesh->NodeVec().Resize(nNodes); //Redimensiona o tamanho do vetor de nos da malha geometrica
  TPZVec<REAL> coord(3,0.0);
  for (long i = 0 ; i < nNodes; i++){
    coord[0] = xCoord[i%(xDiv+1)];
    coord[1] = yCoord[i/(xDiv+1)];

    gmesh->NodeVec()[i].SetCoord(coord); //seta coordenada de um no no vetor de nos da malha
    gmesh->NodeVec()[i].SetNodeId(i);
  }
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
    gmesh->CreateGeoElement(EQuadrilateral, topolQuad, matId, index);//cria elemento quadrilateral
    gmesh->ElementVec()[index];
	}
	
  //CONDICOES DE CONTORNO
  
//  for (int i = 0; i < 2 * xDiv; i++)
//  {
//    topolLine[0] = i % xDiv + (i / xDiv) * yDiv * (xDiv + 1);
//    topolLine[1] = i % xDiv + (i / xDiv) * yDiv * (xDiv + 1) + 1;
//    gmesh->CreateGeoElement(EOned, topolLine, bc0, index);//cria elemento unidimensional
//    gmesh->ElementVec()[index];
//  }
//  
//  for (int i = 0; i < 2 * yDiv; i++)
//  {
//    topolLine[0] = (i % yDiv) * (xDiv + 1) + (i / yDiv ) * xDiv;
//    topolLine[1] = (i % yDiv + 1) * (xDiv + 1) + (i / yDiv ) * xDiv;
//    gmesh->CreateGeoElement(EOned, topolLine, bc0, index);//cria elemento unidimensional
//    gmesh->ElementVec()[index];
//  }
  for (int i = 0; i <xDiv; i++)//baixo
  {
    topolLine[0] = i % xDiv + (i / xDiv) * yDiv * (xDiv + 1);
    topolLine[1] = i % xDiv + (i / xDiv) * yDiv * (xDiv + 1) + 1;
    gmesh->CreateGeoElement(EOned, topolLine, bc0, index);//cria elemento unidimensional
    gmesh->ElementVec()[index];
  }
  for (int i = xDiv; i < 2 * xDiv; i++)//cima
  {
    topolLine[0] = i % xDiv + (i / xDiv) * yDiv * (xDiv + 1);
    topolLine[1] = i % xDiv + (i / xDiv) * yDiv * (xDiv + 1) + 1;
    gmesh->CreateGeoElement(EOned, topolLine, bc1, index);//cria elemento unidimensional
    gmesh->ElementVec()[index];
  }
  
  for (int i = 0; i < yDiv; i++)//esquerda
  {
    topolLine[0] = (i % yDiv) * (xDiv + 1) + (i / yDiv ) * xDiv;
    topolLine[1] = (i % yDiv + 1) * (xDiv + 1) + (i / yDiv ) * xDiv;
    gmesh->CreateGeoElement(EOned, topolLine, bc2, index);//cria elemento unidimensional
    gmesh->ElementVec()[index];
  }
  for (int i = yDiv; i < 2 * yDiv; i++)//direita
  {
    topolLine[0] = (i % yDiv) * (xDiv + 1) + (i / yDiv ) * xDiv;
    topolLine[1] = (i % yDiv + 1) * (xDiv + 1) + (i / yDiv ) * xDiv;
    gmesh->CreateGeoElement(EOned, topolLine, bc3, index);//cria elemento unidimensional
    gmesh->ElementVec()[index];
  }

  
	gmesh->BuildConnectivity(); //constroi a conectividade de vizinhanca da malha
	
	return gmesh;
}


TPZCompMesh *CMesh(TPZGeoMesh *gmesh, int pOrder, STATE (& ur)( TPZVec<REAL>),STATE (& er)( TPZVec<REAL>), REAL freq)
{
  const int dim = 2; //dimensao do problema
  const int dirichlet = 0, neumann = 1; //tipo da condicao de contorno do problema ->default dirichlet na esquerda e na direita
  
  const int matId = 1; //define id para um material(formulacao fraca)
  const int bc0 = -1; //define id para um material(cond contorno baixo)
  const int bc1 = -2; //define id para um material(cond contorno cima)
  const int bc2 = -3; //define id para um material(cond contorno esquerda)
  const int bc3 = -4; //define id para um material(cond contorno direita)
  
  // Criando material
  TPZMatValidacaoHCurlFran1 *material = new TPZMatValidacaoHCurlFran1(matId,freq, ur,er);//criando material que implementa a formulacao fraca do problema de validacao
  
  //AQUIFRAN
  material->SetForcingFunction(forcedRHS);
  ///criar malha computacional
  TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDefaultOrder(pOrder);//seta ordem polimonial de aproximacao
  cmesh->SetDimModel(dim);//seta dimensao do modelo
  // Inserindo material na malha
  cmesh->InsertMaterialObject(material);
		
  ///Inserir condicao de contorno condutores //AQUIFRAN
//  TPZFMatrix<STATE> val1(1,1,0.), val2(1,1,0.);
//  TPZMaterial * BCond0 = material->CreateBC(material, bc0, dirichlet, val1, val2);//cria material que implementa a condicao de contorno de dirichlet
  
  TPZFMatrix<STATE> val1(1,1,0.), val2(1,1,0.);
  //val2(0,0) = -2;
  TPZMaterial * BCond0 = material->CreateBC(material, bc0, dirichlet, val1, val2);//cria material que implementa a condicao de contorno de baixo
  //val2(0,0) = 0.;
  TPZMaterial * BCond1 = material->CreateBC(material, bc1, dirichlet, val1, val2);//cria material que implementa a condicao de contorno de cima
  //val2(0,0) = -2;
  TPZMaterial * BCond2 = material->CreateBC(material, bc2, dirichlet, val1, val2);//cria material que implementa a condicao de contorno da esquerda
  //val2(0,0) = 0.;
  TPZMaterial * BCond3 = material->CreateBC(material, bc3, dirichlet, val1, val2);//cria material que implementa a condicao de contorno da direita

  cmesh->InsertMaterialObject(BCond0);//insere material na malha
  cmesh->InsertMaterialObject(BCond1);//insere material na malha
  cmesh->InsertMaterialObject(BCond2);//insere material na malha
  cmesh->InsertMaterialObject(BCond3);//insere material na malha

  
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



void createRotationMatrix (TPZFNMatrix<9,REAL> &mat, const REAL thetaRotation, const int axis)
{
  mat.Identity();
  switch (axis) {
    case 0://x axis
      mat(1,1) = cos(thetaRotation);
      mat(1,2) = -sin(thetaRotation);
      mat(2,1) = sin(thetaRotation);
      mat(2,2) = cos(thetaRotation);
      break;
    case 1: //y axis
      mat(0,0) = cos(thetaRotation);
      mat(0,2) = -sin(thetaRotation);
      mat(2,0) = sin(thetaRotation);
      mat(2,2) = cos(thetaRotation);
      break;
    case 2://z axis
      mat(0,0) = cos(thetaRotation);
      mat(0,1) = -sin(thetaRotation);
      mat(1,0) = sin(thetaRotation);
      mat(1,1) = cos(thetaRotation);
      break;
    default:
      break;
  }
}



void forcedRHS (const TPZVec<REAL> &xRotated, TPZVec<STATE> &out)
{
  
  //  out[0] = -1. * ( (x[1] - 1.) * (x[1] + 1.) + 2. );
  //  out[1] = -1. * ( (x[0] - 1.) * (x[0] + 1.) + 2. );
  //  out[2] = 0;
  
  
  TPZManVector<REAL,3> x(3,0.);
  x[0] = rotationMatrix(0,0) * xRotated[0] + rotationMatrix(1,0) * xRotated[1] + rotationMatrix(2,0) * xRotated[2];
  x[1] = rotationMatrix(0,1) * xRotated[0] + rotationMatrix(1,1) * xRotated[1] + rotationMatrix(2,1) * xRotated[2];
  x[2] = rotationMatrix(0,2) * xRotated[0] + rotationMatrix(1,2) * xRotated[1] + rotationMatrix(2,2) * xRotated[2];
  TPZManVector<REAL,3> outRotated(3,0.);
  outRotated[0] = -1. * ( (x[1] - 1.) * (x[1] + 1.) + 2. );
  outRotated[1] = -1. * ( (x[0] - 1.) * (x[0] + 1.) + 2. );
  outRotated[2] = 0;
  
  out[0] = rotationMatrix(0,0) * outRotated[0] + rotationMatrix(0,1) * outRotated[1] + rotationMatrix(0,2) * outRotated[2];
  out[1] = rotationMatrix(1,0) * outRotated[0] + rotationMatrix(1,1) * outRotated[1] + rotationMatrix(1,2) * outRotated[2];
  out[2] = rotationMatrix(2,0) * outRotated[0] + rotationMatrix(2,1) * outRotated[1] + rotationMatrix(2,2) * outRotated[2];
  
  
}
