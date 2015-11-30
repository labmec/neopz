/**
 * @file
 * @brief A second contact with complex state-variables using HCurl functions. An electromagnetic wave
 * incides in a dielectric slab with a metal backing. Its
 * reflection coefficient is analised as function of the
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
#include "TPZMatValidacaoHCurlFran2.h"
#include "pzlog.h"
#include "TPZTimer.h"

#include "pzskylstrmatrix.h"
#include "pzstrmatrix.h"
#include "pzstepsolver.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZFrontSym.h"
#include "TPZSkylineNSymStructMatrix.h"


#include "pzgengrid.h"
#include "tpzhierarquicalgrid.h"
#include "pzgeoquad.h"
#include "pzgeopoint.h"
#include "tpzgeoelrefpattern.h"

//------Second Validation of a HCurl Formulation-----------------
enum meshTypeE{ createRectangular=1, createTriangular, createZigZag};

void CreateGMesh(TPZGeoMesh * &gmesh, const int meshType, const REAL hDomain, const REAL wDomain, const REAL x0, const REAL z0, const int xDiv, const int zDiv, int &indexRBC);
/**
 * @brief Creates gmesh corresponding to the simple domain
 * @param hDomain height of the simulation domain
 * @param wDomain width of the simulation domain
 */
TPZGeoMesh *CreateRectangularGMesh(const REAL hDomain, const REAL wDomain, const int xDiv, const int zDiv, int &indexRBC);

void Parametricfunction(const TPZVec<REAL> &par, TPZVec<REAL> &X);

void Parametricfunction2(const TPZVec<REAL> &par, TPZVec<REAL> &X);


TPZGeoMesh *CreateTriangularGMesh(const REAL hDomain, const REAL wDomain, const int xDiv, const int yDiv, int &indexRBC);

void FlipTriangles(TPZGeoMesh *gmesh, int nx);

void FlipTriangles(TPZGeoMesh *gmesh, int nx);

/**
 * @brief Creates gmesh corresponding to the simple domain
 * @param hDomain height of the simulation domain
 * @param wDomain width of the simulation domain
 */
TPZGeoMesh *CreateZigZagGMesh(const REAL hDomain, const REAL wDomain, const int xDiv, const int zDiv, int &indexRBC);

/**
 * @brief Creates cmesh
 * @note All the relevant parameters are arguments of this function (except constants)
 * @param gmesh the geometric mesh
 * @param pOrder polynomial approximation order
 * @param ur relative permeability of the dielectric
 * @param er relative permittivity of the dielectric
 * @param freq frequency of the plane-wave
 */
TPZCompMesh *CMesh(TPZGeoMesh *gmesh, int pOrder, STATE (*ur)( TPZVec<REAL> &),STATE (*er)( TPZVec<REAL> &), REAL freq, REAL theta, REAL e0, REAL lambda, REAL scale);

/**
 * @brief Implements permeability of possibly inhomongeneous media
 * @param x spatial coordinates
 */
inline STATE urSubs(TPZVec<REAL> &x);

/**
 * @brief Implements permittivity of possibly inhomongeneous media
 * @param x spatial coordinates
 */
inline STATE erSubs(TPZVec<REAL> &x);

void ExportCSV(const TPZFMatrix<STATE> matriz, const char* name);

int main(int argc, char *argv[])
{
  HDivPiola = 1;
  TPZTimer timer;
#ifdef LOG4CXX
  InitializePZLOG();
#endif
  
  //PARAMETROS FISICOS DO PROBLEMA
  REAL lambda = 1550*1e-9;
  REAL freq = M_C/lambda;
  REAL theta = 0.;
  REAL e0 = 1.;
  //PARAMETROS DA GEOMETRIA
  REAL hDomain = 5*lambda;
  REAL wDomain = 5*lambda;
  REAL scale = 1.;///(5.*lambda);
  //  REAL hDomain = 1;
  //  REAL wDomain = 1;
  //PARAMETROS DE SIMULACAO
  
  int pOrder = 1; //ordem polinomial de aproximacao
  int dim = 2;
  int xDiv = 1;
  int zDiv = 1;
  
  const int meshType = createTriangular;
  timer.start();
  int indexRBC;
  const REAL x0 = wDomain/2;
  const REAL z0 = hDomain/2;
  
  TPZGeoMesh *gmesh = new TPZGeoMesh();
  CreateGMesh(gmesh, meshType, hDomain, wDomain, x0, z0, xDiv, zDiv, indexRBC);
  
  TPZCompMesh *cmesh = CMesh(gmesh, pOrder, urSubs , erSubs , freq, theta, e0, lambda, scale); //funcao para criar a malha computacional
  bool optimizeBandwidth = false;
  TPZAnalysis an(cmesh,optimizeBandwidth);
  //configuracoes do objeto de analise
  TPZSkylineNSymStructMatrix skylstr(cmesh); //CUIDADO
  skylstr.SetNumThreads(4);
  an.SetStructuralMatrix(skylstr);
  TPZStepSolver<STATE> step;
  step.SetDirect(ELDLt); //caso simetrico
  an.SetSolver(step);
  
  TPZStack<std::string> scalnames, vecnames;
  vecnames.Push("realE");//setando para imprimir campoeletrico
  std::string plotfile= "../ValidacaoHCurlFran2EField.vtk";//arquivo de saida que estara na pasta debug
  an.DefineGraphMesh(dim, scalnames, vecnames, plotfile);//define malha grafica
  int postProcessResolution = 3 ;//define resolucao do pos processamento
  //fim das configuracoes do objeto de analise
  
  int nIteracoes = 1;
  TPZFMatrix<REAL> results(nIteracoes,2);
  //  const REAL k0 = 2*M_PI*freq*sqrt(M_UZERO*M_EZERO);
  //  const REAL L = wDomain ;
  {
    std::ofstream file("../cmesh.txt");
    cmesh->ExpandSolution();
    cmesh->CleanUpUnconnectedNodes();
    cmesh->Print(file);
  }
  
  // Resolvendo o Sistema
  an.Assemble();
  TPZStructMatrix stiff = an.StructMatrix();
  an.Solve();
  const TPZFMatrix<STATE> solucao=cmesh->Solution();//Pegando o vetor de solucao, alphaj
  std::string name("../../sol");
  name.append( std::to_string(meshType) );
  name.append(".csv");
  ExportCSV(solucao,name.c_str());
  TPZFMatrix<STATE> sds(solucao);
  solucao.Print(std::cout);
  
  an.PostProcess(postProcessResolution);//realiza pos processamento*)
  
  for (int i = 0; i < solucao.Rows(); i++) {
    
    sds.Zero();
    sds(i,0) = solucao.GetVal(i,0);
    cmesh->Solution() = sds;
    an.PostProcess(postProcessResolution);//realiza pos processamento*)
  }
  
  timer.stop();
  
  std::cout <<"Tempo de simulacao total = "<<timer.seconds()<<" s\n";
  std::cout << "FINISHED!" << std::endl;
  
  return 0;
}

void CreateGMesh(TPZGeoMesh * &gmesh, const int meshType, const REAL hDomain, const REAL wDomain, const REAL x0, const REAL z0, const int xDiv, const int zDiv, int &indexRBC)
{
  switch (meshType) {
    case createRectangular:
    {
      gmesh = CreateRectangularGMesh(hDomain, wDomain, xDiv, zDiv, indexRBC); //funcao para criar a malha geometrica
      std::ofstream outTxt("../gmeshRectangular.txt"); //define arquivo de saida para impressao da malha no
      gmesh->Print(outTxt);
      std::ofstream out0("../gmeshRectangular.vtk"); //define arquivo de saida para impressao da malha no paraview
      TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out0, true); //imprime a malha no formato vtk
    }
      break;
    case createTriangular:
    {
      gmesh = CreateTriangularGMesh(hDomain, wDomain, xDiv, zDiv, indexRBC); //funcao para criar a malha geometrica
      std::ofstream outTxt("../gmeshTriangular.txt"); //define arquivo de saida para impressao da malha no
      gmesh->Print(outTxt);
      std::ofstream out0("../gmeshTriangular.vtk"); //define arquivo de saida para impressao da malha no paraview
      TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out0, true); //imprime a malha no formato vtk
      
    }
      break;
    case createZigZag:
    {
      gmesh = CreateZigZagGMesh(hDomain, wDomain, xDiv, zDiv, indexRBC);
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

TPZGeoMesh *CreateRectangularGMesh(const REAL hDomain, const REAL wDomain, const int xDiv, const int zDiv, int &indexRBC)
{
  TPZManVector<int,3> nx(3,0);
  TPZManVector<REAL,3> llCoord(3,0.) , ulCoord(3,0.) , urCoord(3,0.) , lrCoord(3,0.);
  llCoord[0] = 0.;
  llCoord[1] = 0.;
  
  ulCoord[0] = 0.;
  ulCoord[1] = hDomain;
  
  urCoord[0] = wDomain;
  urCoord[1] = hDomain;
  
  lrCoord[0] = wDomain;
  lrCoord[1] = 0.;
  
  nx[0]=xDiv;
  nx[1]=zDiv;
  int numl = 1;
  REAL rot = 0.5;
  TPZGenGrid gengrid = TPZGenGrid(nx, llCoord, urCoord,  numl, rot);
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
  gmesh->ResetConnectivities();
  //converte malha xy para xz
  for (int i = 0; i < gmesh->NNodes() ; i ++)
  {
    TPZGeoNode node = gmesh->NodeVec()[i];
    node.SetCoord( 2, node.Coord(1) );
    node.SetCoord( 1, 0. );
    gmesh->NodeVec()[i] = node;
  }
  
  gmesh->BuildConnectivity();
  return gmesh;
  //  TPZGeoMesh * gmesh = new TPZGeoMesh;//Inicializa objeto da classe TPZGeoMesh
  //
  //  int matId = 1; //define id para um material(formulacao fraca)
  //  int bc0 = -1; //define id para um material(cond contorno dirichlet)
  //  int bc1 = -2; //define id para um material(cond contorno mista)
  //  TPZManVector<REAL> xCoord(xDiv+1,0.0);
  //  TPZManVector<REAL> zCoord(zDiv+1,0.0);
  //
  //  int nel = xDiv * zDiv;
  //  int nNodes = (xDiv + 1) * (zDiv + 1);
  //
  //  //coordenadas dos nos da malha sem refinamento
  //  for (int i = 0; i<= xDiv; i++) {
  //    xCoord[i] = 0*(1-(REAL)i/xDiv)+wDomain*((REAL)i/xDiv);
  //  }
  //  for (int i = 0; i<= zDiv; i++) {
  //    zCoord[i] = 0*(1-(REAL)i/zDiv)+hDomain*((REAL)i/zDiv);
  //  }
  //
  //  gmesh->NodeVec().Resize(nNodes); //Redimensiona o tamanho do vetor de nos da malha geometrica
  //  TPZVec<REAL> coord(3,0.0);
  //  for (long i = 0 ; i < nNodes; i++){
  //    coord[0] = xCoord[i%(xDiv+1)];
  //    coord[2] = zCoord[i/(xDiv+1)];
  //
  //    gmesh->NodeVec()[i].SetCoord(coord); //seta coordenada de um no no vetor de nos da malha
  //    gmesh->NodeVec()[i].SetNodeId(i);
  //  }
  //#ifdef DEBUG
  //  for (int i = 0; i< nNodes; i++) {
  //    gmesh->NodeVec()[i].Print();
  //  }
  //#endif
  //  // Criando Elementos
  //  TPZVec <long> topolQuad(4); //vetor que sera inicializado com o indice dos nos de um elemento quadrado bidimensional
  //  TPZVec <long> topolLine(2); //vetor que sera inicializado com o indice dos nos de um elemento unidimensional
  //  long index; //id do elemento que sera preenchido pelo metodo CreateGeoElement
  //  for (long iel = 0; iel < nel; iel++) {
  //    const long ino1 = iel % (xDiv) + (iel / xDiv) * (xDiv + 1);
  //    const long ino2 = iel % (xDiv) + (iel / xDiv + 1) * (xDiv + 1);
  //    const long ino3 = iel % (xDiv) + (iel / xDiv + 1) * (xDiv + 1) + 1;
  //    const long ino4 = iel % (xDiv) + (iel / xDiv) * (xDiv + 1) + 1;
  //    topolQuad[0] = ino1;
  //    topolQuad[1] = ino2;
  //    topolQuad[2] = ino3;
  //    topolQuad[3] = ino4;
  //    //std::cout <<topolQuad[0]<<" "<<topolQuad[1]<<" "<<topolQuad[2]<<" "<<topolQuad[3]<<std::endl;
  //    gmesh->CreateGeoElement(EQuadrilateral, topolQuad, matId, index);//cria elemento quadrilateral
  //    //gmesh->ElementVec()[index];
  //  }
  //
  //  //CONDICOES DE CONTORNO
  //  for (int i = 0; i < 2 * xDiv; i++)
  //  {
  //    topolLine[0] = i % xDiv + (i / xDiv) * zDiv * (xDiv + 1);
  //    topolLine[1] = i % xDiv + (i / xDiv) * zDiv * (xDiv + 1) + 1;
  //    //std::cout <<topolLine[0]<<" "<<topolLine[1]<<std::endl;
  //    gmesh->CreateGeoElement(EOned, topolLine, bc0, index);//cria elemento unidimensional
  //    //gmesh->ElementVec()[index];
  //  }
  //
  //  for (int i = 0; i < 2 * zDiv; i++)
  //  {
  //    topolLine[0] = (i % zDiv) * (xDiv + 1) + (i / zDiv ) * xDiv;
  //    topolLine[1] = (i % zDiv + 1) * (xDiv + 1) + (i / zDiv ) * xDiv;
  //    //std::cout <<topolLine[0]<<" "<<topolLine[1]<<std::endl;
  //    if( !(i / zDiv) )
  //      gmesh->CreateGeoElement(EOned, topolLine, bc0, index);//cria elemento unidimensional
  //    else
  //    {
  //      gmesh->CreateGeoElement(EOned, topolLine, bc1, index);//cria elemento unidimensional
  //      if( i == zDiv )
  //        indexRBC = index;
  //    }
  //    //gmesh->ElementVec()[index];
  //  }
  //
  //
  //  gmesh->BuildConnectivity(); //constroi a conectividade de vizinhanca da malha
  //
  //  return gmesh;
}

TPZGeoMesh *CreateTriangularGMesh(const REAL hDomain, const REAL wDomain, const int xDiv, const int zDiv, int &indexRBC)
{
  TPZGeoMesh * gmesh = new TPZGeoMesh;//Inicializa objeto da classe TPZGeoMesh
  
  int matId = 1; //define id para um material(formulacao fraca)
  int bc0 = -1; //define id para um material(cond contorno dirichlet)
  int bc1 = -2; //define id para um material(cond contorno mista)
  TPZManVector<REAL> xCoord(xDiv+1,0.0);
  TPZManVector<REAL> zCoord(zDiv+1,0.0);
  
  int nNodes = (xDiv + 1) * (zDiv + 1);
  
  //coordenadas dos nos da malha sem refinamento
  for (int i = 0; i<= xDiv; i++) {
    xCoord[i] = 0*(1-(REAL)i/xDiv)+wDomain*((REAL)i/xDiv);
  }
  for (int i = 0; i<= zDiv; i++) {
    zCoord[i] = 0*(1-(REAL)i/zDiv)+hDomain*((REAL)i/zDiv);
  }
  
  gmesh->NodeVec().Resize(nNodes); //Redimensiona o tamanho do vetor de nos da malha geometrica
  TPZVec<REAL> coord(3,0.0);
  for (long i = 0 ; i < nNodes; i++){
    coord[0] = xCoord[i%(xDiv+1)];
    coord[2] = zCoord[i/(xDiv+1)];
    
    gmesh->NodeVec()[i].SetCoord(coord); //seta coordenada de um no no vetor de nos da malha
    gmesh->NodeVec()[i].SetNodeId(i);
  }
  for (int i = 0; i< nNodes; i++) {
    gmesh->NodeVec()[i].Print();
  }
  // Criando Elementos
  TPZVec <long> topolTri(4); //vetor que sera inicializado com o indice dos nos de um elemento quadrado bidimensional
  TPZVec <long> topolLine(2); //vetor que sera inicializado com o indice dos nos de um elemento unidimensional
  long index; //id do elemento que sera preenchido pelo metodo CreateGeoElement
  int nel = xDiv * zDiv * 2;
  for (long iel = 0; iel < nel; iel++) {
    long ino1,ino2,ino3;
    if ((iel+1)%2) {
      ino1 = iel % (xDiv) + (iel / xDiv) * (xDiv + 1);
      ino2 = iel % (xDiv) + (iel / xDiv + 1) * (xDiv + 1) + 1;
      ino3 = iel % (xDiv) + (iel / xDiv + 1) * (xDiv + 1);
    }
    else{
      ino1 = (iel-1) % (xDiv) + ((iel-1) / xDiv) * (xDiv + 1);
      ino2 = (iel-1) % (xDiv) + ((iel-1) / xDiv) * (xDiv + 1) + 1;
      ino3 = (iel-1) % (xDiv) + ((iel-1) / xDiv + 1) * (xDiv + 1) + 1;
      
    }
    topolTri[0] = ino1;
    topolTri[1] = ino2;
    topolTri[2] = ino3;
    gmesh->CreateGeoElement(ETriangle, topolTri, matId, index);//cria elemento triangular
  }
  
  //CONDICOES DE CONTORNO
  for (int i = 0; i < 2 * xDiv; i++)
  {
    topolLine[0] = i % xDiv + (i / xDiv) * zDiv * (xDiv + 1);
    topolLine[1] = i % xDiv + (i / xDiv) * zDiv * (xDiv + 1) + 1;
    //std::cout <<topolLine[0]<<" "<<topolLine[1]<<std::endl;
    gmesh->CreateGeoElement(EOned, topolLine, bc0, index);//cria elemento unidimensional
    //gmesh->ElementVec()[index];
  }
  
  for (int i = 0; i < 2 * zDiv; i++)
  {
    topolLine[0] = (i % zDiv) * (xDiv + 1) + (i / zDiv ) * xDiv;
    topolLine[1] = (i % zDiv + 1) * (xDiv + 1) + (i / zDiv ) * xDiv;
    //std::cout <<topolLine[0]<<" "<<topolLine[1]<<std::endl;
    if( !(i / zDiv) )
      gmesh->CreateGeoElement(EOned, topolLine, bc0, index);//cria elemento unidimensional
    else
    {
      gmesh->CreateGeoElement(EOned, topolLine, bc1, index);//cria elemento unidimensional
      if( i == zDiv )
        indexRBC = index;
    }
  }
  
  
  gmesh->BuildConnectivity(); //constroi a conectividade de vizinhanca da malha
  
  return gmesh;
  
}

TPZGeoMesh *CreateZigZagGMesh(const REAL hDomain, const REAL wDomain, const int xDiv, const int zDiv, int &indexRBC)
{
  TPZManVector<int,3> nx(3,0);
  TPZManVector<REAL,3> llCoord(3,0.) , ulCoord(3,0.) , urCoord(3,0.) , lrCoord(3,0.);
  llCoord[0] = 0.;
  llCoord[1] = 0.;
  
  ulCoord[0] = 0.;
  ulCoord[1] = hDomain;
  
  urCoord[0] = wDomain;
  urCoord[1] = hDomain;
  
  lrCoord[0] = wDomain;
  lrCoord[1] = 0.;
  
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
  
  //converte malha xy para xz
  for (int i = 0; i < gmesh->NNodes() ; i ++)
  {
    TPZGeoNode node = gmesh->NodeVec()[i];
    node.SetCoord( 2, node.Coord(1) );
    node.SetCoord( 1, 0. );
    gmesh->NodeVec()[i] = node;
  }
  
  
  indexRBC = gengrid.ElemId(nx[1]-1,nx[0]-1, numl-1);
  return gmesh;
}

TPZCompMesh *CMesh(TPZGeoMesh *gmesh, int pOrder, STATE (*ur)( TPZVec<REAL> &),STATE (*er)( TPZVec<REAL> &), REAL freq, REAL theta, REAL e0, REAL lambda, REAL scale)
{
  const int dim = 2; //dimensao do problema
  const int matId = 1; //define id para um material(formulacao fraca)
  const int bc0 = -1; //define id para um material(cond contorno dirichlet)
  const int bc1 = -2; //define id para um material(cond contorno mista)
  enum{ dirichlet = 0, neumann, mixed}; //tipo da condicao de contorno do problema
  // Criando material
  TPZMatValidacaoHCurlFran2 *material = new TPZMatValidacaoHCurlFran2(matId,freq, ur,er, theta, scale);//criando material que implementa a formulacao fraca do problema de validacao
  
  ///criar malha computacional
  TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDefaultOrder(pOrder);//seta ordem polimonial de aproximacao
  cmesh->SetDimModel(dim);//seta dimensao do modelo
  // Inserindo material na malha
  cmesh->InsertMaterialObject(material);
		
  ///Inserir condicao de contorno condutores
  TPZFNMatrix<1,STATE> val1(1,1,0.), val2(1,1,0.);
  
  TPZMaterial * BCond0 = material->CreateBC(material, bc0, dirichlet, val1, val2);//cria material que implementa a condicao de contorno de dirichlet
  
  REAL k0 = 2*M_PI*freq*sqrt(M_UZERO*M_EZERO);
  REAL L = 5 * M_C/freq ;
  val1(0,0) = -1.*imaginary*k0*cos(theta);
  val2(0,0) = -2.*imaginary*k0*cos(theta)*e0*exp(imaginary*k0*cos(theta)*L);
  std::cout<<val2<<std::endl;
  TPZMaterial * BCond1 = material->CreateBC(material, bc1, mixed, val1, val2);//cria material que implementa a condicao de contorno mista
  
  cmesh->InsertMaterialObject(BCond0);//insere material na malha
  cmesh->InsertMaterialObject(BCond1);//insere material na malha
  
  cmesh->SetAllCreateFunctionsHDiv();//define espaco de aproximacao
  
  //Cria elementos computacionais que gerenciarao o espaco de aproximacao da malha
  cmesh->AutoBuild();
  
  if (pOrder == 1) {
    //cmesh->CleanUpUnconnectedNodes();
    TPZCreateApproximationSpace::MakeRaviartThomas(*cmesh);
    cmesh->CleanUpUnconnectedNodes();
  }
  //cmesh->AutoBuild();
  return cmesh;
}

inline STATE urSubs(TPZVec<REAL> &x)
{
  return ( 2.-imaginary*0.1 );
}

inline STATE erSubs(TPZVec<REAL> &x)
{
  return ( 4.+(2.-imaginary*0.1)*(1.-x[0]/(5*1550*1e-9))*(1.-x[0]/(5*1550*1e-9)) );
}
void ExportCSV(const TPZFMatrix<STATE> matriz, const char* name)
{
  std::ofstream file;
  if(name!=NULL) file.open(name);
  
  for (int i = 0 ; i < matriz.Rows(); i++) {
    for(int j = 0 ; j < matriz.Cols() ; j++){
      file<<matriz.GetVal(i, j);
      if (j != matriz.Cols() - 1) {
        file<<" , ";
      }
    }
    file<<std::endl;
  }
  file.close();
}

