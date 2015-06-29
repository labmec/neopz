/**
 * @file
 * @brief A First contact with complex state-variables. An electromagnetic wave
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
#include "TPZMatComplexExample2D.h"
#include "TPZTimer.h"

//------Plane-wave reflection by a metal-backed dielectric slab-----------------



/**
 * @brief Creates gmesh with nel unidimensional elements with elsize length
 * @param nel number of elements
 * @param elsize element length
 * @param indexRBC index of the element correspondent to the right boundary of the domain
 */
TPZGeoMesh *CreateGMesh(long nelX, long nelZ, REAL L, REAL H, TPZVec <long> &indexRBC);

/**
 * @brief Creates cmesh
 * @note All the relevant parameters are arguments of this function (except constants)
 * @param gmesh the geometric mesh
 * @param pOrder polynomial approximation order
 * @param L length of the dielectric slab
 * @param ur relative permeability of the dielectric
 * @param er relative permittivity of the dielectric
 * @theta angle of incidence at the slab
 * @lambda wavelength of the plane-wave
 * @e0 magnitude of the electric field
 */
TPZCompMesh *CMesh(TPZGeoMesh *gmesh, int pOrder, REAL L, STATE (& ur)( TPZVec<REAL>, REAL),STATE (& er)( TPZVec<REAL>,REAL), REAL theta, REAL lambda,REAL e0);

/**
 * @brief implements the relative permeability of the material as a function of position and wavelength
 * @param x spatial coordinates
 * @param lambda wavelength of the wave in the dielectric
 */
STATE urMat(TPZVec<REAL> x, REAL lambda);

/**
 * @brief implements the relative permittivity of the material as a function of position and wavelength
 * @param x spatial coordinates
 * @param lambda wavelength of the wave in the dielectric
 */
STATE erMat(TPZVec<REAL> x, REAL lambda);

/**
 * @brief Export a TPZFMatrix instance as a mathematica file with the matrix and a listplot.nb
 * @param matrix desired matrix to be exported
 * @param matrix_name desired name of the variable ibn the nb
 * @param nel number of elements (to identify the plot)
 * @param plot_name PlotLabel->"plot_name"
 */
void PrintMathematica(TPZFMatrix<REAL> matrix, const char *matrix_name, int nelX, int nelZ, const char *plot_name);


using namespace std;



int main(int argc, char *argv[])
{
  TPZTimer timer;
  
  //PARAMETROS FISICOS DO PROBLEMA
  REAL lambda = 1550*1e-9;//comprimento de onda
  REAL e0=1.;//amplitude do campo
  REAL w=2.*M_PI*M_C/lambda;
  REAL kZero=w*sqrt(M_UZERO*M_EZERO);
  
  REAL L = 5*lambda; //largura do slab com inicio na origem zero
  REAL H = 5*lambda; //domínio será um quadrado em x-z
  int nelX = 3; //numero de divisoes em x
  int nelZ = 2; //numero de divisoes em x
  int pOrder = 1; //ordem polinomial de aproximacao
  
  
//	REAL elsize = L/nel; //tamanho do lado de cada elemento
  TPZVec <long> indexRBC(nelZ,0);
  
  int nIteracoes=100;
  std::cout<<"Numero de iteracoes = "<<nIteracoes<<std::endl;
  std::cout<<"Numero de elementos = "<<nelX*nelZ<<std::endl;
  std::cout<<"Ordem polinomial de aproximacao = "<<pOrder<<std::endl;
  
  
  timer.start();
  
	TPZGeoMesh *gmesh = CreateGMesh(nelX, nelZ, L, H, indexRBC); //funcao para criar a malha geometrica
	
	std::ofstream out("gmesh.vtk"); //define arquivo de saida para impressao da malha no paraview
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true); //imprime a malha no formato vtk
  gmesh->Print();
	
  TPZFMatrix<REAL> results(nIteracoes,2);
  
  REAL theta = 0.;
  TPZCompMesh *cmesh = CMesh(gmesh, pOrder, L, urMat,erMat, theta, lambda,e0); //funcao para criar a malha computacional
  bool optimizeBandwidth = true; //impede a renumeracao das equacoes do problema(para obter o mesmo resultado do Oden)
  TPZAnalysis an(cmesh,optimizeBandwidth);
  
  for(int i=0;i<nIteracoes;i++)
  {

    // Resolvendo o Sistema
    an.Run();//assembla a matriz de rigidez (e o vetor de carga) global e inverte o sistema de equacoes
    
    //Calcula coeficiente de reflexao para z = ?
    TPZManVector<REAL,3> qsi(3,0.);
    TPZManVector<STATE,3> sol(3,0.);
    int var = 0;
    TPZCompEl *cel = cmesh->Element(indexRBC[0]);
    if (!cel){
      DebugStop();
    }
    cel->Solution(qsi, var, sol);
    STATE eZApprox = sol[0]+imaginary*sol[1];
    
    REAL R=abs((eZApprox-e0*exp(imaginary*kZero*L*cos(theta)))/(e0*exp(-1.*imaginary*kZero*L*cos(theta))));
    results(i,0)=theta*180/M_PI;
    results(i,1)=R;
    
    //Atualiza theta para a proxima iteracao
    theta=((i+1.)/(nIteracoes-1))*M_PI/2.;
    TPZMaterial *mat = cmesh->FindMaterial(1);
    TPZMatComplexExample2D *mpmat = dynamic_cast<TPZMatComplexExample2D *>(mat);
    mpmat->SetTheta(theta);
    
    //Atualiza condicoes de contorno
    const STATE val1=-1.*imaginary*kZero*cos(theta);
    const STATE val2=-2.*imaginary*kZero*cos(theta)*e0*exp(imaginary*kZero*L*cos(theta));
    
    TPZBndCond *mpmatbc = dynamic_cast<TPZBndCond *>(cmesh->FindMaterial(-2));
    mpmatbc->Val1()(0,0) = val1;
    mpmatbc->Val2()(0,0) = val2;
  }

  timer.stop();
  
  PrintMathematica(results, "resultadoPZ",nelX,nelZ,"Plot do modulo de reflexao em funcao do angulo de incidencia");
  //fazendo pos processamento para paraview
//  TPZStack<string> scalnames, vecnames;
//  scalnames.Push("State");//setando para imprimir u
//  string plotfile= "ComplexExample2DSol.vtk";//arquivo de saida que estara na pasta debug
//  int dim=1;
//  an.DefineGraphMesh(dim, scalnames, vecnames, plotfile);//define malha grafica
//  int postProcessResolution = 0;//define resolucao do pos processamento
//  an.PostProcess(postProcessResolution);//realiza pos processamento
  

  std::cout <<"Tempo de simulacao total = "<<timer.seconds()<<" s\n";
	std::cout << "FINISHED!" << std::endl;
	
	return 0;
}



TPZGeoMesh *CreateGMesh(long nelX, long nelZ, REAL L, REAL H, TPZVec <long> &indexRBC)
{
	TPZGeoMesh * gmesh = new TPZGeoMesh;//Inicializa objeto da classe TPZGeoMesh
	
	long nnodes = (nelX+1)*(nelZ+1); //numero de nos do problema
	gmesh->NodeVec().Resize(nnodes); //Redimensiona o tamanho do vetor de nos da malha geometrica
	int mat1d = 1; //define id para um material(formulacao fraca)
  int bc0 = -1; //define id para um material(cond contorno esq)
  int bc1 = -2; //define id para um material(cond contorno dir)
  int bc2 = -3; //define id para um material (cond contorno top/bottom)
	
  REAL incX = L / nelX, //largura dos elementos
  incY = H / nelZ; //altura dos elementos
  // Colocando nos na malha
  for (long i = 0 ; i < nnodes; i++){
    long iX = i / (nelZ + 1);
    long iY = i % (nelZ + 1);
    TPZVec <REAL> coord(3,0.);
    coord[0] = iX * incX;
    coord[1] = iY * incY;
    
    gmesh->NodeVec()[i].SetCoord(coord); //seta coordenada de um no no vetor de nos da malha
    gmesh->NodeVec()[i].SetNodeId(i);
  }
	
	// Criando Elementos
  TPZVec <long> topolQuad(4); //vetor que sera inicializado com o indice dos nos de um elemento quadrado bidimensional
	TPZVec <long> topolLine(2); //vetor que sera inicializado com o indice dos nos de um elemento unidimensional
//  TPZVec <long> TopolPoint(1); //vetor que sera inicializado com o indice do no de um elemento zero-dimensional
	long index; //id do elemento que sera preenchido pelo metodo CreateGeoElement
  long nel = nelX*nelZ;
	for (long iel = 0; iel < nel; iel++) {
		const long ino1 = iel % nelZ + (iel / nelZ) * ( nelZ + 1 ) ;//LL
		const long ino2 = iel % nelZ + (iel / nelZ) * ( nelZ + 1 ) + 1;//UL
    const long ino3 = iel % nelZ + (iel / nelZ + 1) * ( nelZ + 1 ) + 1;//UR
    const long ino4 = iel % nelZ + (iel / nelZ + 1) * ( nelZ + 1 );//LR
		topolQuad[0] = ino1;
		topolQuad[1] = ino2;
    topolQuad[2] = ino3;
    topolQuad[3] = ino4;
    gmesh->CreateGeoElement(EQuadrilateral, topolQuad, mat1d, index);//cria elemento quadrilateral
    gmesh->ElementVec()[index];
	}
	
	// Cond Contorno esquerda
  for (long iBC = 0; iBC < nelZ; iBC++){
    const long ino1 = iBC;
    const long ino2 = iBC + 1;
    topolLine[0] = ino1;
    topolLine[1] = ino2;
    gmesh->CreateGeoElement(EOned, topolLine, bc0, index);//cria elemento unidimensional
    gmesh->ElementVec()[index];
  }

	// Cond Contorno Direita
  for (long iBC = 0; iBC < nelZ; iBC++) {
    const long ino1 = (nelZ + 1) * nelX + iBC;
    const long ino2 = (nelZ + 1) * nelX + iBC + 1;
    topolLine[0] = ino1;
    topolLine[1] = ino2;
    gmesh->CreateGeoBlendElement(EOned, topolLine, bc1, index);//cria elemento unidimensional
    gmesh->ElementVec()[index];
    indexRBC[iBC] = index;
  }
  
  // Cond Contorno Fundo
  for (long iBC = 0; iBC < nelX; iBC++) {
    const long ino1 = iBC * (nelZ + 1);
    const long ino2 = (iBC + 1) * (nelZ + 1);
    topolLine[0] = ino1;
    topolLine[1] = ino2;
    gmesh->CreateGeoBlendElement(EOned, topolLine, bc2, index);//cria elemento unidimensional
    gmesh->ElementVec()[index];
  }
  
  // Cond Contorno Topo
  for (long iBC = 0; iBC < nelX; iBC++) {
    const long ino1 = (iBC + 1) * (nelZ + 1) - 1;
    const long ino2 = (iBC + 2) * (nelZ + 1) - 1;
    topolLine[0] = ino1;
    topolLine[1] = ino2;
    gmesh->CreateGeoBlendElement(EOned, topolLine, bc2, index);//cria elemento unidimensional
    gmesh->ElementVec()[index];
  }
  
	gmesh->BuildConnectivity(); //constroi a conectividade de vizinhanca da malha
	
	return gmesh;
}

TPZCompMesh *CMesh(TPZGeoMesh *gmesh, int pOrder, REAL L, STATE (& ur)( TPZVec<REAL>, REAL),STATE (& er)( TPZVec<REAL>,REAL), REAL theta, REAL lambda,REAL e0)
{
	const int dim = 2; //dimensao do problema
	const int matId = 1, bc0 = -1, bc1 = -2; //MESMOS ids da malha geometrica
	const int dirichlet = 0, neumann = 1, mixed = 2; //tipo da condicao de contorno do problema ->default dirichlet na esquerda e na direita 
	
	// Criando material
  //TPZMatComplexExample2D(int id, REAL l, REAL t, REAL eO,STATE (* ur)( TPZVec<REAL>),STATE (* er)( TPZVec<REAL>));
	TPZMatComplexExample2D *material = new TPZMatComplexExample2D(matId,lambda,theta,e0,ur,er);//criando material que implementa a formulacao fraca do problema modelo
	
	///criar malha computacional
	TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDefaultOrder(pOrder);//seta ordem polimonial de aproximacao
	cmesh->SetDimModel(dim);//seta dimensao do modelo
	
	// Inserindo material na malha
	cmesh->InsertMaterialObject(material);
		
	///Inserir condicao de contorno esquerda
	TPZFMatrix<STATE> val1(1,1,0.), val2(1,1,0.);
	TPZMaterial * BCond0 = material->CreateBC(material, bc0, dirichlet, val1, val2);//cria material que implementa a condicao de contorno da esquerda
	
	// Condicao de contorno da direita

  REAL w=2.*M_PI*M_C/lambda;
  REAL kZero=w*sqrt(M_UZERO*M_EZERO);
  
  val1(0,0)=-1.*imaginary*kZero*cos(theta);
  val2(0,0)=-2.*imaginary*kZero*cos(theta)*e0*exp(imaginary*kZero*L*cos(theta));
	TPZMaterial * BCond1 = material->CreateBC(material, bc1, mixed, val1, val2);//cria material que implementa a condicao de contorno da direita
	
	cmesh->InsertMaterialObject(BCond0);//insere material na malha
	cmesh->InsertMaterialObject(BCond1);//insere material na malha
  
  cmesh->SetAllCreateFunctionsHDiv();//define espaco de aproximacao
	
	//Cria elementos computacionais que gerenciarao o espaco de aproximacao da malha
	cmesh->AutoBuild();
	
	return cmesh;
	
}

//implementar condicao de contorno

STATE urMat(TPZVec<REAL> x, REAL l)
{
  STATE val(0.);
  val=2.-imaginary*0.1;
  return val;
}

STATE erMat(TPZVec<REAL> x, REAL l)
{
  STATE val(0.);
  REAL pos=x[0];
  REAL L=5*l;
  val=4.+(2.-imaginary*0.1)*(1.-pos/L)*(1.-pos/L);
  return val;
}

void PrintMathematica(TPZFMatrix<REAL> matriz, const char *nome_matriz,int nelX,int nelZ, const char *titulo_plot)
{
  std::string str="";
  stringstream ss;
  ss<<nome_matriz;
  str.append(ss.str());
  str.append(".nb");
  char * name = new char[str.size() + 1];
  std::copy(str.begin(), str.end(), name);
  name[str.size()] = '\0'; // don't forget the terminating 0
  
  // don't forget to free the string after finished using it
  ofstream file (name, ios::out);
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
      file<<matriz(i,j);
      
    }
    file<<"}";
  }
  file<<"};"<<std::endl;
  
  file<<"ListPlot["<<nome_matriz<<",Joined->True,ImageSize->700,PlotMarkers->{Automatic,6},PlotLegends->{\""<<nelX<<"x"<<nelZ<<" elementos\"},PlotLabel->\""<<titulo_plot<<"\"]"<<std::endl;

}