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


/**
 * @brief Creates gmesh corresponding to the simple domain
 * @param hDomain height of the simulation domain
 * @param wDomain width of the simulation domain
 */
TPZGeoMesh *CreateRectangularGMesh(REAL hDomain, REAL wDomain, int xDiv, int zDiv, int &indexRBC);

void Parametricfunction(const TPZVec<REAL> &par, TPZVec<REAL> &X);

void Parametricfunction2(const TPZVec<REAL> &par, TPZVec<REAL> &X);

TPZGeoMesh *CreateTriangularGMesh(REAL hDomain, REAL wDomain, REAL x0, REAL y0, int xDiv, int yDiv, int &indexRBC);

/**
 * @brief Creates gmesh corresponding to the simple domain
 * @param hDomain height of the simulation domain
 * @param wDomain width of the simulation domain
 */
TPZGeoMesh *CreateZigZagGMesh(REAL hDomain, REAL wDomain, int xDiv, int zDiv, int &indexRBC);

/**
 * @brief Creates cmesh
 * @note All the relevant parameters are arguments of this function (except constants)
 * @param gmesh the geometric mesh
 * @param pOrder polynomial approximation order
 * @param ur relative permeability of the dielectric
 * @param er relative permittivity of the dielectric
 * @param freq frequency of the plane-wave
 */
TPZCompMesh *CMesh(TPZGeoMesh *gmesh, int pOrder, STATE (& ur)( TPZVec<REAL>),STATE (& er)( TPZVec<REAL>), REAL freq, REAL theta, REAL e0, REAL lambda);

/**
 * @brief Implements permeability of possibly inhomongeneous media
 * @param x spatial coordinates
 */
inline STATE urSubs(TPZVec<REAL> x);

/**
 * @brief Implements permittivity of possibly inhomongeneous media
 * @param x spatial coordinates
 */
inline STATE erSubs(TPZVec<REAL> x);

void PrintMathematica(TPZFMatrix<STATE> matrix, const char *matrix_name, int nel, const char *plot_name);

REAL kScale = 1.;
int __attribute__((flatten)) main(int argc, char *argv[])
{
  TPZTimer timer;
#ifdef LOG4CXX
  InitializePZLOG();
#endif
  
  //PARAMETROS FISICOS DO PROBLEMA
  REAL lambda = kScale * 1550*1e-9;
  REAL freq = M_C/lambda;
  REAL theta = 0.;
  REAL e0 = 1.;
  //PARAMETROS DA GEOMETRIA
  REAL hDomain = 5*lambda;
  REAL wDomain = 5*lambda;
//  REAL hDomain = 1;
//  REAL wDomain = 1;
  //PARAMETROS DE SIMULACAO
  
  int pOrder = 1; //ordem polinomial de aproximacao
  int dim = 2;
  int xDiv = 1;
  int zDiv = 1;
  enum meshType{ createRectangular, createTriangular, createZigZag};
  const int mesh = createRectangular;
  timer.start();
  int indexRBC;
  const REAL x0 = wDomain/2;
  const REAL z0 = hDomain/2;

  TPZGeoMesh *gmesh = new TPZGeoMesh();
  
  switch (mesh) {
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
      gmesh = CreateTriangularGMesh(hDomain, wDomain, x0, z0, xDiv, zDiv, indexRBC); //funcao para criar a malha geometrica
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

#ifdef DEBUG 

  gmesh->Print();
#endif
  
  TPZCompMesh *cmesh = CMesh(gmesh, pOrder, urSubs , erSubs , freq, theta, e0,lambda); //funcao para criar a malha computacional
  bool optimizeBandwidth = false;
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
  std::string plotfile= "../ValidacaoHCurlFran2EField.vtk";//arquivo de saida que estara na pasta debug
  an.DefineGraphMesh(dim, scalnames, vecnames, plotfile);//define malha grafica
  int postProcessResolution = 4 ;//define resolucao do pos processamento
  //fim das configuracoes do objeto de analise
  
  int nIteracoes = 1;
  TPZFMatrix<REAL> results(nIteracoes,2);
  const REAL k0 = 2*M_PI*freq*sqrt(M_UZERO*M_EZERO);
  const REAL L = wDomain ;
  {
    std::ofstream file("../cmesh.txt");
    cmesh->Print(file);
  }
  for(int i=0;i<nIteracoes;i++)
  {
    
    // Resolvendo o Sistema
//    an.Assemble();
//    TPZMatrixSolver<STATE> &matrix = an.Solver();
//    TPZFMatrix<STATE> fLoad(16, 16, 0.), result(16,16,0.);
//    
//    fLoad.Identity();
//    matrix.Matrix()->Multiply(fLoad,result);
//    result.Print("rhs",std::cout,EMathematicaInput);
    an.Run();//assembla a matriz de rigidez (e o vetor de carga) global e inverte o sistema de equacoes

      //fazendo pos processamento para paraview (theta = 0)
#ifdef DEBUG
    std::cout<<"realizando pos processamento"<<std::endl;
#endif
    an.PostProcess(postProcessResolution);//realiza pos processamento*)
#ifdef DEBUG
    std::cout<<"finalizando pos processamento"<<std::endl;
#endif
    
//      TPZFMatrix<STATE> solucao=cmesh->Solution();//Pegando o vetor de solucao, alphaj
//      TPZFNMatrix<200,STATE> sds(solucao);
//      for (int i = 0; i < solucao.Rows(); i++) {
//        sds.Zero();
//        sds(i,0) = 1.;
//        cmesh->Solution() = sds;
//        an.PostProcess(postProcessResolution);//realiza pos processamento*)
//      }
//    sds.Zero();
//    sds(1,0) = 1.;
//    sds(0,0) = 1.;
//    cmesh->Solution() = sds;
//    an.PostProcess(postProcessResolution);//realiza pos processamento*)
//    
//    sds.Zero();
//    sds(3,0) = 1.;
//    sds(2,0) = 1.;
//    cmesh->Solution() = sds;
//    an.PostProcess(postProcessResolution);
//    sds.Zero();
//    sds(5,0) = 1.;
//    sds(4,0) = 1.;
//    cmesh->Solution() = sds;
//    an.PostProcess(postProcessResolution);
    //    TPZFMatrix<STATE> solucao=cmesh->Solution();//Pegando o vetor de solucao, alphaj
    //    STATE eZApprox = solucao(solucao.Rows()-1,0);

//    TPZCompEl *cel = cmesh->Element(indexRBC);
//    if (!cel){
//      DebugStop();
//    }
//    TPZManVector<REAL,1> qsi(1,0.);//DUVIDA
//    TPZManVector<STATE,3> sol(3,0.);
//    int var = 0;
//    cel->Solution(qsi, var, sol);
//    STATE eZApprox = 1. * sol[0]+ 1. * imaginary*sol[1];
//    std::cout<<std::endl<<"theta "<<theta<<std::endl<<"sol[0]"<<sol[0]<<std::endl<<"sol[1]"<<sol[1]<<std::endl<<"|eZApprox|"<<std::abs(eZApprox)<<std::endl;
//    
//    STATE gamma =(eZApprox-e0*exp(imaginary*k0*L*cos(theta)))/(e0*exp(-1.*imaginary*k0*L*cos(theta)));
//    std::cout<<"complex reflection coefficient "<<gamma<<std::endl;
//    REAL R = std::abs(gamma);
//    results(i,0)=theta*180/M_PI;
//    results(i,1)=R;
    
//    //avanca para proxima iteracao
//    theta=((i+1.)/(nIteracoes-1))*M_PI/2.;
//    
//    TPZMatValidacaoHCurlFran2 *mpmat = dynamic_cast<TPZMatValidacaoHCurlFran2 *>(cmesh->FindMaterial(1));
//    mpmat->SetTheta(theta);
//    const STATE val1=1.*imaginary*k0*cos(theta);
//    const STATE val2=2.*imaginary*k0*cos(theta)*e0*exp(imaginary*k0*cos(theta)*L);
//    
//    TPZBndCond *mpmatbc = dynamic_cast<TPZBndCond *>(cmesh->FindMaterial(-2));
//    mpmatbc->Val1()(0,0) = val1;
//    mpmatbc->Val2()(0,0) = val2;
  }
  
  timer.stop();
  
//    TPZFMatrix<STATE> solucao=cmesh->Solution();//Pegando o vetor de solucao, alphaj
//    TPZFNMatrix<200,STATE> sds(solucao);
//    for (int i = 0; i < sds.Rows(); i++) {
//      sds.Zero();
//      sds(i,0) = 1.;
//      cmesh->Solution() = sds;
//      an.PostProcess(postProcessResolution);//realiza pos processamento*)
//    }
  
//  PrintMathematica(results, "resultadoPZ",xDiv,"Plot do modulo de reflexao em funcao do angulo de incidencia");
//  an.Run();
//  
//  
//  //an.Rhs().Print("FelPZ",std::cout,EMathematicaInput);
//  timer.stop();
  
  
  std::cout <<"Tempo de simulacao total = "<<timer.seconds()<<" s\n";
  std::cout << "FINISHED!" << std::endl;
  
  return 0;
}



TPZGeoMesh *CreateRectangularGMesh(REAL hDomain, REAL wDomain, int xDiv, int zDiv, int &indexRBC)
{
  TPZGeoMesh * gmesh = new TPZGeoMesh;//Inicializa objeto da classe TPZGeoMesh
  
  int matId = 1; //define id para um material(formulacao fraca)
  int bc0 = -1; //define id para um material(cond contorno dirichlet)
  int bc1 = -2; //define id para um material(cond contorno mista)
  TPZManVector<REAL> xCoord(xDiv+1,0.0);
  TPZManVector<REAL> zCoord(zDiv+1,0.0);
  
  int nel = xDiv * zDiv;
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
    //gmesh->ElementVec()[index];
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
    //gmesh->ElementVec()[index];
  }
  
  
  gmesh->BuildConnectivity(); //constroi a conectividade de vizinhanca da malha
  
  return gmesh;
}

TPZGeoMesh *CreateTriangularGMesh(REAL hDomain, REAL wDomain, REAL x0, REAL y0,int xDiv, int yDiv, int &indexRBC)
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
  
//  teste invariancia da solucao em relacao a reoordenacao dos nos do elemento
//  TPZAdmChunkVector<TPZGeoEl *> elVec =  gmesh->ElementVec();
//  TPZManVector<int,3> nodeIndices;
//  int hasChanged = 0;
//  for (int i = 0; i < elVec.NElements() ; i++) {
//    
//    elVec[i]->GetNodeIndices(nodeIndices);
//    if (nodeIndices.size() != 3) {
//      continue;
//    }
//    if (hasChanged == 0) {
//      hasChanged = 1;
//      continue;
//    }
//    TPZManVector<int,3> nodeIndicesCp(nodeIndices);
//    for (int j = 0; j < nodeIndices.size(); j++) {
//      nodeIndices[j] = nodeIndicesCp[(j+1) % 3];
//      elVec[i]->SetNodeIndex(j, nodeIndices[j]);
//    }
//    
//  }
//  
//  gmesh->ResetConnectivities();
//  gmesh->BuildConnectivity();
  
  
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
  
  return gmesh;

}

TPZGeoMesh *CreateZigZagGMesh(REAL hDomain, REAL wDomain, int xDiv, int zDiv, int &indexRBC)
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

TPZCompMesh *CMesh(TPZGeoMesh *gmesh, int pOrder, STATE (& ur)( TPZVec<REAL>),STATE (& er)( TPZVec<REAL>), REAL freq, REAL theta, REAL e0, REAL lambda)
{
  const int dim = 2; //dimensao do problema
  const int matId = 1; //define id para um material(formulacao fraca)
  const int bc0 = -1; //define id para um material(cond contorno dirichlet)
  const int bc1 = -2; //define id para um material(cond contorno mista)
  const int dirichlet = 0, mixed = 2; //tipo da condicao de contorno do problema
  // Criando material
  TPZMatValidacaoHCurlFran2 *material = new TPZMatValidacaoHCurlFran2(matId,freq, ur,er, theta);//criando material que implementa a formulacao fraca do problema de validacao
  
  ///criar malha computacional
  TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDefaultOrder(pOrder);//seta ordem polimonial de aproximacao
  cmesh->SetDimModel(dim);//seta dimensao do modelo
  // Inserindo material na malha
  cmesh->InsertMaterialObject(material);
		
  ///Inserir condicao de contorno condutores
  TPZFNMatrix<1,STATE> val1(1,1,0.), val2(1,1,0.);
  
  val1(0,0) = 1.e22/kScale;
  TPZMaterial * BCond0 = material->CreateBC(material, bc0, mixed, val1, val2);//cria material que implementa a condicao de contorno de dirichlet
  
  REAL k0 = 2*M_PI*freq*sqrt(M_UZERO*M_EZERO);
  REAL L = 5 * M_C/freq ;
  val1(0,0) = 1.*imaginary*k0*cos(theta);
  val2(0,0) = 2.*imaginary*k0*cos(theta)*e0*exp(imaginary*k0*cos(theta)*L);
  TPZMaterial * BCond1 = material->CreateBC(material, bc1, mixed, val1, val2);//cria material que implementa a condicao de contorno mista
  
  cmesh->InsertMaterialObject(BCond0);//insere material na malha
  cmesh->InsertMaterialObject(BCond1);//insere material na malha
  
  cmesh->SetAllCreateFunctionsHDiv();//define espaco de aproximacao
  
  //Cria elementos computacionais que gerenciarao o espaco de aproximacao da malha
  cmesh->AutoBuild();
  
  return cmesh;
}

inline STATE urSubs(TPZVec<REAL> x)
{
  return ( 2.-imaginary*0.1 );
}

inline STATE erSubs(TPZVec<REAL> x)
{
  return ( 4.+(2.-imaginary*0.1)*(1.-x[0]/(5*kScale*1550*1e-9))*(1.-x[0]/(5*kScale*1550*1e-9)) );
}

void PrintMathematica(TPZFMatrix<STATE> matriz, const char *nome_matriz,int nel,const char *titulo_plot)
{
  std::string str="";
  std::stringstream ss;
  ss<<nome_matriz;
  str.append(ss.str());
  str.append(".nb");
  char * name = new char[str.size() + 1];
  std::copy(str.begin(), str.end(), name);
  name[str.size()] = '\0'; // don't forget the terminating 0
  
  // don't forget to free the string after finished using it
  std::ofstream file (name, std::ios::out);
  delete[] name;
  if(!file.is_open())
    DebugStop();
  file<<nome_matriz<<"=";
  file<<"{";
  for(int i=0;i<matriz.Rows();i++)
  {
    if(i>0)
    {
      file<<", ";
    }
    file<<"{";
    for(int j=0;j<matriz.Cols();j++)
    {
      if(j>0)
      {
        file<<", ";
      }
      file<<"{"<<std::real(matriz(i,j))<<"+I*"<<std::imag(matriz(i,j))<<"}";
      
    }
    file<<"}";
  }
  file<<"};"<<std::endl;
  
  file<<"ListPlot["<<nome_matriz<<",Joined->True,ImageSize->700,PlotMarkers->{Automatic,6},PlotLegends->{\""<<nel<<" elementos\"},PlotLabel->\""<<titulo_plot<<"\"]"<<std::endl;
  
}