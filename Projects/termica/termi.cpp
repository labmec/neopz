//Classes utilitárias
#include "pzvec.h"
#include "pzstack.h"
#include "TPZTimer.h"

//Classes Geométricas
#include "pzgmesh.h"
#include "pzgeoel.h"

//Classes Computacionais
#include "pzcmesh.h"
#include "pzcompel.h"

//Materiais
#include "pzelasmat.h"
#include "pzmat2dlin.h"
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

//Análise
#include "pzanalysis.h"

//Reordenamento de equações grafos
#include "pzmetis.h"
#include "tpznodesetcompute.h"

//Pós-processamento
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

//Método para inserção de materiais
void InicializarMaterial(TPZCompMesh &cmesh);

//Método para a leitura da malha
void LerMalha(char *arquivo,TPZGeoMesh &mesh);

//Método para a criação de um elemento
void UmElemento(TPZGeoMesh &malha);

//Inserção de uma condição de contorno dada por uma função
void forcingfunction(TPZVec<REAL> &ponto, TPZVec<REAL> &force);

int main() {
  //Arquivo de impressão de log (não essencial)
  std::ofstream file("graphtest.txt");
   
  //malha geometrica
  TPZGeoMesh *malha = new TPZGeoMesh;
  malha->SetName("Malha geometrica Curso Floripa");
   
  //Cria e insere um elemento na malha geométrica criada acima
  UmElemento(*malha);
  
  // Os elementos da malha original podem ser refinados
  TPZVec<int> csub(0);
  TPZManVector<TPZGeoEl *> pv(4);  //vetor auxiliar que receberá os subelementos 
  int n1=1,level=0;   
  cout << "\nDividir ate nivel ? ";
  int resp = 3;
  cin >> resp;      
  
  int nelc = malha->ElementVec().NElements(); //número de elementos da malha
  int el;  //iterador

  TPZGeoEl *cpel;  //ponteiro auxiliar para um elemento geométrico
  for(el=0;el<malha->ElementVec().NElements();el++) { //loop sobre todos os elementos da malha
    cpel = malha->ElementVec()[el];    // obtém o ponteiro para o elemento da malha de 
                                       // posição igual ao do iterador atual
    if(cpel && cpel->Level() < resp){  // verifica se o ponteiro não é nulo e se
                                       // se o nível de refinamento é menor que o nível desejado
      cpel->Divide(pv);                // Executa a divisão do elemento
    }
  }
  
  //solicita ao usuário uma ordem p de interpolacao padrão para a malha
  int ord = 4;
  cout << "Entre ordem 1,2,3,4,5 : ";
  cin >> ord;
  
  //Define a ordem p de criação de todo elemento computacional como sendo ord
  TPZCompEl::gOrder = ord;
  
  // Construção da malha computacional. Aqui a malha terá apenas uma referência
  // para a malha geométrica...
  TPZCompMesh *malhacomp = new TPZCompMesh(malha);
  
  //Pode-se dar um nome a malha para efeito de impressão de resultados...
  malhacomp->SetName("Malha Computacional :  Connects e Elementos");
  
  // Define-se o problema a resolver em cada parte do domínio / definem-se os materiais
  InicializarMaterial(*malhacomp);
  
  //Define-se o espaço de aproximação
  malhacomp->AutoBuild();
  
  //Ajusta o espaço levando em conta as condições de contorno
  malhacomp->AdjustBoundaryElements();
  
  //Inicializa-se a estrutura de dados
  malhacomp->InitializeBlock();
  
   
  //Verificação e melhoria da banda da matriz de rigidez (Opcional)
  TPZFMatrix before, after;
  malhacomp->ComputeFillIn(100,before);

  // Cria-se um objeto análise. O único parâmetro é a malha computacional
  // com toda a sua estrutura de dados já preparada. O PZ tentará otimizar
  // a banda aqui.
  TPZAnalysis an(malhacomp);
  
  // cálculo da banda após o ajuste
  malhacomp->ComputeFillIn(100,after);
  
  //Pós processamento da matriz de rigidez antes e depois da otimização da banda 
  VisualMatrix(before,"before.dx");
  VisualMatrix(after,"after.dx");
   
  //Escolha do padrão de armazenamento. O parâmetro de entrada é a malha computacional
  //TPZSkylineStructMatrix strmat(malhacomp);
  TPZParSkylineStructMatrix strmat(malhacomp);       //skyline em paralelo (multthread)
  //TPZParFrontStructMatrix<TPZFrontSym> strmat(malhacomp);
  
  //Define-se o padrão de armazenamento para a análise
  an.SetStructuralMatrix(strmat);
  
  //Criação e definição do solver 
  TPZStepSolver step;
  step.SetDirect(ECholesky);   
  
  //Define-se o solver para a análise
  an.SetSolver(step);
  /*
  TPZMatrixSolver *precond = an.BuildPreconditioner(TPZAnalysis::EElement,false);
  step.SetCG(100,*precond,1.e-10,0);
  an.SetSolver(step);
  delete precond;
  */
  //Processamento
  cout << "Number of equations " << malhacomp->NEquations() << endl;
  an.Run();

  // Posprocessamento
  TPZVec<char *> scalnames(1);
  scalnames[0] = "state";    //nome das variáveis que se quer pós-processar
  TPZVec<char *> vecnames(0);
  char plotfile[] =  "termica.dx"; //nome do arquivo de resultados
  //Define-se para a análise as variáveis a pós-processar
  an.DefineGraphMesh(2, scalnames, vecnames, plotfile);
  //Executa os cálculos para geração dos resultados de pós-processamento
  an.PostProcess(2);
  //   an.DefineGraphMesh(2, scalnames, vecnames, pltfile);
  //   an.PostProcess(2);
  
  //limpar a memória alocada dinamicamente.
  delete malhacomp;
  delete malha;
  return 0;
}

void UmElemento(TPZGeoMesh &malha) {
  //Para efeito de teste será criado um elemento quadrilateral de
  //comprimento 1, com as seguintes coordenadas (x,y,z)
  // - canto inferior esquerdo: (0,0,0);
  // - canto inferior direito : (0,1,0);
  // - canto superior direito : (1,1,0);
  // - canto superior esquerdo: (0,1,0);  
  double coordstore[4][3] = {{0.,0.,0.},{1.,0.,0.},{1.,1.,0.},{0.,1.,0.}};

  // criar os quatro nós geométricos
  int i,j;                              //iteradores
  TPZVec<REAL> coord(3,0.);             //vetor auxiliar para armazenar uma coordenada
  
  for(i=0; i<4; i++) {                  //loop sobre o número de nós da malha
    
    // inicializar as coordenadas do no em um vetor do tipo pz
    for (j=0; j<3; j++) coord[j] = coordstore[i][j]; 
    // identificar um espaço no vetor da malha geométrica 
    // onde podemos armazenar o objeto nó a criar
    int nodeindex = malha.NodeVec ().AllocateNewElement ();

    // criar um nó geométrico e inserí-lo na posição 
    // alocada no vetor de nós da malha geométrica
    malha.NodeVec ()[i].Initialize (i,coord,malha);
  
  }

  // Criação de um elemento geométrico
  // inicializar os índices dos nós do elemento
  TPZVec<int> indices(4);
  for(i=0; i<4; i++) indices[i] = i; //loop sobre o número de nós do elemento
  //no caso só há quatro nós e eles foram criados na ordem correta para o
  //elemento em questão. A ordem dos nós deve seguir um padrão pré-estabelecido
  
  // O próprio construtor vai inserir o elemento na malha
  // os parâmetros de criação do elemento são os seguintes:
  // 1) Tipo geométrico do elemento
  // - EPoint            element 0D - type point        -  associated index 0
  // - EOned             element 1D - type oned         -  associated index 1
  // - ETriangle         element 2D - type triangle     -  associated index 2
  // - EQuadrilateral    element 2D - type quad         -  associated index 3
  // - ETetraedro        element 3D - type tetraedro    -  associated index 4
  // - EPiramide         element 3D - type piramide     -  associated index 5
  // - EPrisma           element 3D - type prisma       -  associated index 6
  // - ECube             element 3D - type cube         -  associated index 7
  // 2) Vetor de conectividades dos elementos (o número de nós deve ser 
  //    compatível com o número de nós do tipo do elemento
  // 3) Índice do material que será associado ao elemento
  // 4) Variável onde será retornado o índice do elemento criado no vetor
  //    de elementos da malha geométrica
  // TPZGeoEl *gel = new TPZGeoElQ2d(0,indices,1,malha);  //forma antiga.
  int index;
  TPZGeoEl *gel = malha.CreateGeoElement(EQuadrilateral,indices,1,index);

  //Gerar as estruturas de dados de conectividade e vizinhança
  malha.BuildConnectivity ();

  // Identificar onde serão inseridas condições de contorno.
  // Uma condição de contorno é aplicada a um lado (parte do contorno) 
  // de um elemento. Este objeto irá inserir-se automaticamente na malha
  // Os parâmetros são os seguintes:
  // 1) elemento onde será aplicada a condição de contorno
  // 2) lado do elemento onde será inserida a condição de contorno
  // 3) identificador da condição de contorno
  // 4) referência para a malha geométrica.
  TPZGeoElBC(gel,4,-4,malha);
}


void forcingfunction(TPZVec<REAL> &p, TPZVec<REAL> &f) {
  REAL x0=0.5,y0=0.2,r=0.005,eps=0.0001;

  REAL dist2 = (p[0]-x0)*(p[0]-x0)+(p[1]-y0)*(p[1]-y0);
  REAL A = exp(-dist2/eps);
  REAL Bx2 = 4.*(p[0]-x0)*(p[0]-x0)/(eps*eps);
  REAL By2 = 4.*(p[1]-y0)*(p[1]-y0)/(eps*eps);
  REAL Bx = 2.*(p[0]-x0);
  REAL By = 2.*(p[1]-y0);
  if(dist2 < r*r) {
    f[0] = A*Bx2-2.*A/eps+A*By2-2.*A/eps;
  } else f[0] = 0.;
}

// Exemplo de leitura de malha a partir de um arquivo
// os parâmetros utilizados são o nome do arquivo e uma
// referência para a malha geométrica que receberá os dados
void LerMalha(char *nome, TPZGeoMesh &grid) {
  //Criação de um objeto arquivo de entrada de dados;
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
  //O parâmetro básico para a criação de um material é o 
  // identificador. Esse deve corresponder aos índices utilizados
  // na criação dos elementos geométricos.
  // No restante, o material define a equação diferencial a ser
  // resolvida...
  TPZMat2dLin *meumat = new TPZMat2dLin(1);
  
  //Cada material tem parâmetros de inicialização próprios, assim
  //deve-se consultar a documentação para verificar como definir os
  //parâmetros. No caso em questão o material requer três matrizes
  //e uma função de cálculo também é fornecida
  TPZFMatrix xk(1,1,1.),xc(1,2,0.),xf(1,1,1.);
  meumat->SetMaterial (xk,xc,xf);
  meumat->SetForcingFunction(forcingfunction);
  
  //Após a criação do material este dever ser inserido na estrutura
  //de dados da malha computacional
  cmesh.InsertMaterialObject(meumat);

  // inserir as condições de contorno
  // Uma condição de contorno pode ser dada por duas matrizes
  // relacionadas a rigidez e ao vetor de carga. As dimensões dessas
  // matrizes dependem da dimensão e do número de variáveis de estado
  // do problema. O tipo de condição de contorno é definido por um 
  // identificador, onde:
  //  - 0 = Dirichlet
  //  - 1 = Neumann
  //  - 2 = Mista
  TPZFMatrix val1(1,1,0.),val2(1,1,0.);
  
  //Os parâmetros necessários à criação de uma condição de contorno são:
  // 1) Identificador da condição de contorno (lembre-se que uma BC é como um material)
  // 2) O Tipo da BC : Dirichlet, Neumann ou Mista 
  // 3) Os valores da BC
  TPZMaterial *bnd = meumat->CreateBC (-4,0,val1,val2);
  
  //Da mesma forma que para os materiais, após sua criação é necessário
  //a sua inserção na estrutura de dados da malha computacional 
  cmesh.InsertMaterialObject(bnd);
  
  //criação e inserção de outras BC's
  bnd = meumat->CreateBC (-3,0,val1,val2);
  cmesh.InsertMaterialObject(bnd);
  bnd = meumat->CreateBC (-2,0,val1,val2);
  cmesh.InsertMaterialObject(bnd);
  bnd = meumat->CreateBC (-1,0,val1,val2);
  cmesh.InsertMaterialObject(bnd);
}

