//Classes utilitï¿½ias
#include "pzvec.h"
#include "pzstack.h"
#include "TPZTimer.h"

//Classes Geomï¿½ricas
#include "pzgmesh.h"
#include "pzgeoel.h"

//Classes Computacionais
#include "pzcmesh.h"
#include "pzcompel.h"

//Materiais
#include "pzelasmat.h"
#include "pzmat2dlin.h"
#include "pzpoisson3d.h"
#include "pzbndcond.h"

//Matrizes de armazenamento e estruturais
#include "pzfmatrix.h"
//#include "pzskylmat.h"
#include "pzysmp.h"
#include "pzskylstrmatrix.h"
#include "TPZParSkylineStructMatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZFrontSym.h"

//Solver
#include "pzstepsolver.h"

//Anï¿½ise
#include "pzanalysis.h"

//Reordenamento de equaï¿½es grafos
#include "pzmetis.h"
#include "tpznodesetcompute.h"

//Pï¿½-processamento
#include "pzvisualmatrix.h"

//Bibliotecas std, math etc
#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "fast.h"
//#include "pzelct2d.h"
//template<class T>
//class TPZVec;
//#define NOTDEBUG

using namespace std;
//Método para inserção de materiais
void InicializarMaterial(TPZCompMesh &cmesh);

//Mï¿½odo para a leitura da malha
void LerMalha(char *arquivo,TPZGeoMesh &mesh);

//Mï¿½odo para a criaï¿½o de um elemento
void UmElemento(TPZGeoMesh &malha);
void UmElemento3D(TPZGeoMesh &malha);

//Inserï¿½o de uma condiï¿½o de contorno dada por uma funï¿½o
void forcingfunction(TPZVec<REAL> &ponto, TPZVec<REAL> &force);

int main() {
  //Arquivo de impressï¿½ de log (nï¿½ essencial)
  std::ofstream file("graphtest.txt");
   
  //malha geometrica
  TPZGeoMesh *malha = new TPZGeoMesh;
  malha->SetName("Malha geometrica Curso Floripa");
   
  //Cria e insere um elemento na malha geomï¿½rica criada acima
  UmElemento3D(*malha);
  
  // Os elementos da malha original podem ser refinados
  TPZVec<int> csub(0);
  TPZManVector<TPZGeoEl *> pv(4);  //vetor auxiliar que receberï¿½os subelementos 
//  int n1=1,level=0;   
  cout << "\nDividir ate nivel ? ";
  int resp = 1;
//  cin >> resp;      
  
//  int nelc = malha->ElementVec().NElements(); //nmero de elementos da malha
  int el;  //iterador

  TPZGeoEl *cpel;  //ponteiro auxiliar para um elemento geomï¿½rico
  for(el=0;el<malha->ElementVec().NElements();el++) { //loop sobre todos os elementos da malha
    cpel = malha->ElementVec()[el];    // obtï¿½ o ponteiro para o elemento da malha de 
                                       // posiï¿½o igual ao do iterador atual
    if(cpel && cpel->Level() < resp){  // verifica se o ponteiro nï¿½ ï¿½nulo e se
                                       // se o nï¿½el de refinamento ï¿½menor que o nï¿½el desejado
      cpel->Divide(pv);                // Executa a divisï¿½ do elemento
    }
  }
  
  //solicita ao usuï¿½io uma ordem p de interpolacao padrï¿½ para a malha
  int ord = 1;
  cout << "Entre ordem 1,2,3,4,5 : ";
//  cin >> ord;
  
  //Define a ordem p de criaï¿½o de todo elemento computacional como sendo ord
  TPZCompEl::gOrder = ord;
  
  // Construï¿½o da malha computacional. Aqui a malha terï¿½apenas uma referï¿½cia
  // para a malha geomï¿½rica...
  TPZCompMesh *malhacomp = new TPZCompMesh(malha);
  
  //Pode-se dar um nome a malha para efeito de impressï¿½ de resultados...
  malhacomp->SetName("Malha Computacional :  Connects e Elementos");
  
  // Define-se o problema a resolver em cada parte do domï¿½io / definem-se os materiais
  InicializarMaterial(*malhacomp);
  
//  string filename("CursoAltoDes.in");
//  string filename("matriz1089.in");
//  string filename("matriz12.in");
//  string filename("matriz16000.in");
//  string filename("matriz260000.in");
  TPZFMatrix rhs,sol;
  //Define-se o espaï¿½ de aproximaï¿½o
  malhacomp->AutoBuild();
  
  //Ajusta o espaï¿½ levando em conta as condiï¿½es de contorno
  malhacomp->AdjustBoundaryElements();
  
  //Inicializa-se a estrutura de dados
  malhacomp->InitializeBlock();
  
   
  //Verificaï¿½o e melhoria da banda da matriz de rigidez (Opcional)
  TPZFMatrix before, after;
  malhacomp->ComputeFillIn(100,before);

  // Cria-se um objeto anï¿½ise. O nico parï¿½etro ï¿½a malha computacional
  // com toda a sua estrutura de dados jï¿½preparada. O PZ tentarï¿½otimizar
  // a banda aqui.
  TPZAnalysis an(malhacomp);
  
  // cï¿½culo da banda apï¿½ o ajuste
  malhacomp->ComputeFillIn(100,after);
  
  //Pï¿½ processamento da matriz de rigidez antes e depois da otimizaï¿½o da banda 
  VisualMatrix(before,"before.dx");
  VisualMatrix(after,"after.dx");
   
  //Escolha do padrï¿½ de armazenamento. O parï¿½etro de entrada ï¿½a malha computacional
  //TPZSkylineStructMatrix strmat(malhacomp);
  TPZSkylineStructMatrix strmat(malhacomp);       //skyline em paralelo (multthread)
  //TPZParFrontStructMatrix<TPZFrontSym> strmat(malhacomp);
  
  //Define-se o padrï¿½ de armazenamento para a anï¿½ise
  an.SetStructuralMatrix(strmat);
  
  //Criaï¿½o e definiï¿½o do solver 
  TPZStepSolver step;
//  step.SetDirect(ECholesky);   
  //Define-se o solver para a anï¿½ise
  an.SetSolver(step);
  
//  an.Assemble();
  //TPZMatrixSolver *precond = an.BuildPreconditioner(TPZAnalysis::EBlockJacobi,false);
  TPZStepSolver jac;
  jac.SetJacobi(1,0.,0);
  jac.ShareMatrix(step);
  step.SetCG(10000,jac,1.e-10,0);
  an.SetSolver(step);
//  delete precond;
  
  //Processamento
  cout << "Number of equations " << malhacomp->NEquations() << endl;
  an.Run();

  std::ofstream out("output.txt");
  an.Solution().Print("Solution obtained",out);

  // Posprocessamento
  TPZVec<char *> scalnames(1);
  scalnames[0] = "state";    //nome das variï¿½eis que se quer pï¿½-processar
  TPZVec<char *> vecnames(0);
  char plotfile[] =  "termica.dx"; //nome do arquivo de resultados
  //Define-se para a anï¿½ise as variï¿½eis a pï¿½-processar
  an.DefineGraphMesh(2, scalnames, vecnames, plotfile);
  //Executa os cï¿½culos para geraï¿½o dos resultados de pï¿½-processamento
  an.PostProcess(2);
  //   an.DefineGraphMesh(2, scalnames, vecnames, pltfile);
  //   an.PostProcess(2);
  
  //limpar a memï¿½ia alocada dinamicamente.
  delete malhacomp;
  delete malha;
  return 0;
}

void UmElemento(TPZGeoMesh &malha) {
  //Para efeito de teste serï¿½criado um elemento quadrilateral de
  //comprimento 1, com as seguintes coordenadas (x,y,z)
  // - canto inferior esquerdo: (0,0,0);
  // - canto inferior direito : (0,1,0);
  // - canto superior direito : (1,1,0);
  // - canto superior esquerdo: (0,1,0);  
  double coordstore[4][3] = {{0.,0.,0.},{1.,0.,0.},{1.,1.,0.},{0.,1.,0.}};

  // criar os quatro nï¿½ geomï¿½ricos
  int i,j;                              //iteradores
  TPZVec<REAL> coord(3,0.);             //vetor auxiliar para armazenar uma coordenada
  malha.NodeVec().Resize(4);
  
  for(i=0; i<4; i++) {                  //loop sobre o nmero de nï¿½ da malha
    
    // inicializar as coordenadas do no em um vetor do tipo pz
    for (j=0; j<3; j++) coord[j] = coordstore[i][j]; 
    // identificar um espaï¿½ no vetor da malha geomï¿½rica 
    // onde podemos armazenar o objeto nï¿½a criar
//    int nodeindex = malha.NodeVec ().AllocateNewElement ();

    // criar um nï¿½geomï¿½rico e inserï¿½lo na posiï¿½o 
    // alocada no vetor de nï¿½ da malha geomï¿½rica
    malha.NodeVec ()[i].Initialize (i,coord,malha);
  
  }

  // Criaï¿½o de um elemento geomï¿½rico
  // inicializar os ï¿½dices dos nï¿½ do elemento
  TPZVec<int> indices(4);
  for(i=0; i<4; i++) indices[i] = i; //loop sobre o nmero de nï¿½ do elemento
  //no caso sï¿½hï¿½quatro nï¿½ e eles foram criados na ordem correta para o
  //elemento em questï¿½. A ordem dos nï¿½ deve seguir um padrï¿½ prï¿½estabelecido
  
  // O prï¿½rio construtor vai inserir o elemento na malha
  // os parï¿½etros de criaï¿½o do elemento sï¿½ os seguintes:
  // 1) Tipo geomï¿½rico do elemento
  // - EPoint            element 0D - type point        -  associated index 0
  // - EOned             element 1D - type oned         -  associated index 1
  // - ETriangle         element 2D - type triangle     -  associated index 2
  // - EQuadrilateral    element 2D - type quad         -  associated index 3
  // - ETetraedro        element 3D - type tetraedro    -  associated index 4
  // - EPiramide         element 3D - type piramide     -  associated index 5
  // - EPrisma           element 3D - type prisma       -  associated index 6
  // - ECube             element 3D - type cube         -  associated index 7
  // 2) Vetor de conectividades dos elementos (o nmero de nï¿½ deve ser 
  //    compatï¿½el com o nmero de nï¿½ do tipo do elemento
  // 3) ï¿½dice do material que serï¿½associado ao elemento
  // 4) Variï¿½el onde serï¿½retornado o ï¿½dice do elemento criado no vetor
  //    de elementos da malha geomï¿½rica
  // TPZGeoEl *gel = new TPZGeoElQ2d(0,indices,1,malha);  //forma antiga.
  int index;
  TPZGeoEl *gel = malha.CreateGeoElement(EQuadrilateral,indices,1,index);

  //Gerar as estruturas de dados de conectividade e vizinhanï¿½
  malha.BuildConnectivity ();

  // Identificar onde serï¿½ inseridas condiï¿½es de contorno.
  // Uma condiï¿½o de contorno ï¿½aplicada a um lado (parte do contorno) 
  // de um elemento. Este objeto irï¿½inserir-se automaticamente na malha
  // Os parï¿½etros sï¿½ os seguintes:
  // 1) elemento onde serï¿½aplicada a condiï¿½o de contorno
  // 2) lado do elemento onde serï¿½inserida a condiï¿½o de contorno
  // 3) identificador da condiï¿½o de contorno
  // 4) referï¿½cia para a malha geomï¿½rica.
  TPZGeoElBC(gel,4,-4,malha);
}

void UmElemento3D(TPZGeoMesh &malha) {
  //Para efeito de teste serï¿½criado um elemento quadrilateral de
  //comprimento 1, com as seguintes coordenadas (x,y,z)
  // - canto inferior esquerdo: (0,0,0);
  // - canto inferior direito : (0,1,0);
  // - canto superior direito : (1,1,0);
  // - canto superior esquerdo: (0,1,0);  
  double coordstore[8][3] = {{0.,0.,0.},{1.,0.,0.},{1.,1.,0.},{0.,1.,0.},
	{0.,0.,1.},{1.,0.,1.},{1.,1.,1.},{0.,1.,1.}};

  // criar os quatro nï¿½ geomï¿½ricos
  int i,j;                              //iteradores
  TPZVec<REAL> coord(3,0.);             //vetor auxiliar para armazenar uma coordenada
  malha.NodeVec().Resize(8);
  
  for(i=0; i<8; i++) {                  //loop sobre o nmero de nï¿½ da malha
    
    // inicializar as coordenadas do no em um vetor do tipo pz
    for (j=0; j<3; j++) coord[j] = coordstore[i][j]; 
    // identificar um espaï¿½ no vetor da malha geomï¿½rica 
    // onde podemos armazenar o objeto nï¿½a criar
//    int nodeindex = malha.NodeVec ().AllocateNewElement ();

    // criar um nï¿½geomï¿½rico e inserï¿½lo na posiï¿½o 
    // alocada no vetor de nï¿½ da malha geomï¿½rica
    malha.NodeVec ()[i].Initialize (i,coord,malha);
  
  }

  // Criaï¿½o de um elemento geomï¿½rico
  // inicializar os ï¿½dices dos nï¿½ do elemento
  TPZVec<int> indices(8);
  for(i=0; i<8; i++) indices[i] = i; //loop sobre o nmero de nï¿½ do elemento
  //no caso sï¿½hï¿½quatro nï¿½ e eles foram criados na ordem correta para o
  //elemento em questï¿½. A ordem dos nï¿½ deve seguir um padrï¿½ prï¿½estabelecido
  
  // O prï¿½rio construtor vai inserir o elemento na malha
  // os parï¿½etros de criaï¿½o do elemento sï¿½ os seguintes:
  // 1) Tipo geomï¿½rico do elemento
  // - EPoint            element 0D - type point        -  associated index 0
  // - EOned             element 1D - type oned         -  associated index 1
  // - ETriangle         element 2D - type triangle     -  associated index 2
  // - EQuadrilateral    element 2D - type quad         -  associated index 3
  // - ETetraedro        element 3D - type tetraedro    -  associated index 4
  // - EPiramide         element 3D - type piramide     -  associated index 5
  // - EPrisma           element 3D - type prisma       -  associated index 6
  // - ECube             element 3D - type cube         -  associated index 7
  // 2) Vetor de conectividades dos elementos (o nmero de nï¿½ deve ser 
  //    compatï¿½el com o nmero de nï¿½ do tipo do elemento
  // 3) ï¿½dice do material que serï¿½associado ao elemento
  // 4) Variï¿½el onde serï¿½retornado o ï¿½dice do elemento criado no vetor
  //    de elementos da malha geomï¿½rica
  // TPZGeoEl *gel = new TPZGeoElQ2d(0,indices,1,malha);  //forma antiga.
  int index;
  TPZGeoEl *gel = malha.CreateGeoElement(ECube,indices,1,index);

  //Gerar as estruturas de dados de conectividade e vizinhanï¿½
  malha.BuildConnectivity ();

  // Identificar onde serï¿½ inseridas condiï¿½es de contorno.
  // Uma condiï¿½o de contorno ï¿½aplicada a um lado (parte do contorno) 
  // de um elemento. Este objeto irï¿½inserir-se automaticamente na malha
  // Os parï¿½etros sï¿½ os seguintes:
  // 1) elemento onde serï¿½aplicada a condiï¿½o de contorno
  // 2) lado do elemento onde serï¿½inserida a condiï¿½o de contorno
  // 3) identificador da condiï¿½o de contorno
  // 4) referï¿½cia para a malha geomï¿½rica.
  TPZGeoElBC(gel,21,-4,malha);
}


void forcingfunction(TPZVec<REAL> &p, TPZVec<REAL> &f) {
//  REAL x0=0.5,y0=0.2,r=0.005,eps=0.0001;
  REAL x0=0.5,y0=0.4,r=0.2,eps=0.1;

  REAL dist2 = (p[0]-x0)*(p[0]-x0)+(p[1]-y0)*(p[1]-y0);
  REAL A = exp(-dist2/eps);
  REAL Bx2 = 4.*(p[0]-x0)*(p[0]-x0)/(eps*eps);
  REAL By2 = 4.*(p[1]-y0)*(p[1]-y0)/(eps*eps);
//  REAL Bx = 2.*(p[0]-x0);
//  REAL By = 2.*(p[1]-y0);
  if(dist2 < r*r) {
    f[0] = A*Bx2-2.*A/eps+A*By2-2.*A/eps;
  } else f[0] = 0.;
}

// Exemplo de leitura de malha a partir de um arquivo
// os parï¿½etros utilizados sï¿½ o nome do arquivo e uma
// referï¿½cia para a malha geomï¿½rica que receberï¿½os dados
void LerMalha(char *nome, TPZGeoMesh &grid) {
  //Criaï¿½o de um objeto arquivo de entrada de dados;
  ifstream infile(nome);

	int linestoskip;
	char buf[256];
	infile >> linestoskip;
	int i,j;
	for(i=0; i<linestoskip;i++) infile.getline(buf,255);
	infile.getline (buf,255);
	infile.getline (buf,255);
	int ntri,npoin,nbouf,nquad,nsidif;
	infile >> ntri >> npoin >> nbouf >> nquad >> nsidif;
	infile.getline (buf,255);
	infile.getline(buf,255);

	grid.NodeVec ().Resize(npoin+1);
	TPZVec<int> nodeindices(4);
	int mat, elid;
	for(i=0;i<nquad;i++) {
		infile >> elid;
		for(j=0; j<4;j++) infile >> nodeindices[j];
		infile >> mat;
                int index;
                grid.CreateGeoElement(EQuadrilateral,nodeindices,mat,index);
//		new TPZGeoElQ2d(elid,nodeindices,mat,grid);
	}
	infile.getline(buf,255);
	infile.getline(buf,255);

	int nodeid,dum;
	char c;
	TPZVec<REAL> coord(3,0.);
	for(i=0; i<npoin; i++) {
		infile >> nodeid >> coord[0] >> coord[1] >> c >> dum;
		grid.NodeVec ()[nodeid].Initialize (nodeid,coord,grid);
	}
	infile.getline (buf,255);
	infile.getline (buf,255);

	TPZVec<int> sideid(2,0);
	for(i=0; i<nbouf; i++) {
		infile >> sideid[0] >> sideid[1] >> elid >> dum >> mat;
		TPZGeoEl *el = grid.ElementVec ()[elid-1];
		int side = el->WhichSide (sideid);
		TPZGeoElBC(el,side,-mat,grid);
	}
	grid.BuildConnectivity();

	return;
}

void InicializarMaterial(TPZCompMesh &cmesh) {
  //O parï¿½etro bï¿½ico para a criaï¿½o de um material ï¿½o 
  // identificador. Esse deve corresponder aos ï¿½dices utilizados
  // na criaï¿½o dos elementos geomï¿½ricos.
  // No restante, o material define a equaï¿½o diferencial a ser
  // resolvida...

  int dim = cmesh.Reference()->ElementVec()[0]->Dimension();
  TPZMat2dLin *meumat = new TPZMat2dLin(2);
  TPZMatPoisson3d *poismat = new TPZMatPoisson3d(1,dim);
  TPZVec<REAL> convdir(3,0.);
  poismat->SetParameters(1.,0.,convdir);
  poismat->SetInternalFlux(1.);
  //poismat->SetForcingFunction(forcingfunction);
  
  //Cada material tem parï¿½etros de inicializaï¿½o prï¿½rios, assim
  //deve-se consultar a documentaï¿½o para verificar como definir os
  //parï¿½etros. No caso em questï¿½ o material requer trï¿½ matrizes
  //e uma funï¿½o de cï¿½culo tambï¿½ ï¿½fornecida
  TPZFMatrix xk(1,1,1.),xc(1,2,0.),xf(1,1,1.);
  meumat->SetMaterial (xk,xc,xf);
//  meumat->SetForcingFunction(forcingfunction);
  
  //Apï¿½ a criaï¿½o do material este dever ser inserido na estrutura
  //de dados da malha computacional
  cmesh.InsertMaterialObject(meumat);
  cmesh.InsertMaterialObject(poismat);

//  TPZMaterial *atual = meumat;
  TPZMaterial *atual = poismat;

  // inserir as condiï¿½es de contorno
  // Uma condiï¿½o de contorno pode ser dada por duas matrizes
  // relacionadas a rigidez e ao vetor de carga. As dimensï¿½s dessas
  // matrizes dependem da dimensï¿½ e do nmero de variï¿½eis de estado
  // do problema. O tipo de condiï¿½o de contorno ï¿½definido por um 
  // identificador, onde:
  //  - 0 = Dirichlet
  //  - 1 = Neumann
  //  - 2 = Mista
  TPZFMatrix val1(1,1,0.),val2(1,1,0.);
  
  //Os parï¿½etros necessï¿½ios ï¿½criaï¿½o de uma condiï¿½o de contorno sï¿½:
  // 1) Identificador da condiï¿½o de contorno (lembre-se que uma BC ï¿½como um material)
  // 2) O Tipo da BC : Dirichlet, Neumann ou Mista 
  // 3) Os valores da BC
  TPZMaterial *bnd = atual->CreateBC (-4,0,val1,val2);
  
  //Da mesma forma que para os materiais, apï¿½ sua criaï¿½o ï¿½necessï¿½io
  //a sua inserï¿½o na estrutura de dados da malha computacional 
  cmesh.InsertMaterialObject(bnd);
  
  //criaï¿½o e inserï¿½o de outras BC's
  bnd = atual->CreateBC (-3,0,val1,val2);
  cmesh.InsertMaterialObject(bnd);
  bnd = atual->CreateBC (-2,0,val1,val2);
  cmesh.InsertMaterialObject(bnd);
  bnd = atual->CreateBC (-1,0,val1,val2);
  cmesh.InsertMaterialObject(bnd);
}

