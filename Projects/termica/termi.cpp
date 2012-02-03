//Classes utilit�ias
#include "pzvec.h"
#include "pzstack.h"
#include "TPZTimer.h"

//Classes Geom�ricas
#include "pzgmesh.h"
#include "pzgeoel.h"

//Classes Computacionais
#include "pzcmesh.h"
#include "pzcompel.h"

//Materiais
#include "pzelasmat.h"
#include "pzmat2dlin.h"
#include "pzpoisson3d.h"
#include "pzelast3d.h"
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

//An�ise
#include "pzanalysis.h"

//Reordenamento de equa�es grafos
#include "pzmetis.h"
#include "tpznodesetcompute.h"

//P�-processamento
#include "pzvisualmatrix.h"

//Bibliotecas std, math etc
#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "fast.h"

#include "pzgengrid.h"

#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif
using namespace std;
//M�todo para inser��o de materiais
void InicializarMaterial(TPZCompMesh &cmesh);

//M�odo para a leitura da malha
void LerMalha(char *arquivo,TPZGeoMesh &mesh);

//M�odo para a cria�o de um elemento
void UmElemento(TPZGeoMesh &malha);
void UmElemento3D(TPZGeoMesh &malha);

//Inser�o de uma condi�o de contorno dada por uma fun�o
void forcingfunction(TPZVec<REAL> &ponto, TPZVec<REAL> &force);

int main()
{
//    TPZGenGrid(TPZVec<int> &nx, TPZVec<REAL> &x0, TPZVec<REAL> &x1, int numl = 1, REAL rot = 0.5);
    TPZVec<int> nx(2,64);
    TPZVec<REAL> x0(2,0.),x1(2,1.);
    TPZGenGrid grid(nx,x0,x1);
    TPZGeoMesh gmesh;
    grid.Read(gmesh);
    grid.SetBC(&gmesh,0,-1);
    grid.SetBC(&gmesh,1,-1);
    grid.SetBC(&gmesh,2,-1);
    grid.SetBC(&gmesh,3,-1);
    TPZCompMesh cmesh(&gmesh);
    TPZMatPoisson3d *poiss = new TPZMatPoisson3d(1,2);
    TPZAutoPointer<TPZMaterial> autopoiss(poiss);
    TPZFMatrix val1(1,1,0.),val2(1,1,0.);
    cmesh.InsertMaterialObject(autopoiss);
    TPZBndCond *bc = new TPZBndCond(autopoiss, -1, 0, val1, val2);
    TPZAutoPointer<TPZMaterial> bcauto(bc);
    cmesh.InsertMaterialObject(bcauto);
    cmesh.SetDefaultOrder(2);
    cmesh.AutoBuild();
    TPZAnalysis an(&cmesh);
    TPZSkylineStructMatrix strskyl(&cmesh);
    strskyl.SetNumThreads(2);
    an.SetStructuralMatrix(strskyl);
    TPZStepSolver solve;
    solve.SetDirect(ECholesky);
    an.SetSolver(solve);
    std::cout << "before running\n";
    std::cout << "neq " << cmesh.NEquations() << std::endl;
#ifdef USING_BOOST
    boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::local_time();
#endif
    an.Run();
#ifdef USING_BOOST
    boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();
    std::cout << "t1 " << t1 << " t2 " << t2 << " elapse " << t2-t1 << "finished\n";
#endif
    TPZFMatrix vismat;
    cmesh.ComputeFillIn(100,vismat);
    VisualMatrixVTK(vismat,"matrixstruct.vtk");
    return 0;
    
}

int main44(){
	int n = 34, m = 42, p = 25;
	std::ofstream matrices("matrices.txt");
	TPZFMatrix A(m,n);
	double val = 0.;
	for(int i = 0; i < m; i++) for(int j = 0; j < n; j++){
		val = rand()/100000000.;
		A(i,j) = val;
	}
	TPZFMatrix x(n,p);
	for(int i = 0; i < n; i++) for(int j = 0; j < p; j++) x(i,j) = random()/100000000.;
	A.Print("A=", matrices, EMathematicaInput);
	x.Print("x=", matrices, EMathematicaInput);
	TPZFMatrix B;
	A.ConstMultiply(x,B);
	B.Print("B=", matrices, EMathematicaInput);
	return 0;
}

int main22() {
  //Arquivo de impress� de log (n� essencial)
  std::ofstream file("graphtest.txt");

  //malha geometrica
  TPZGeoMesh *malha = new TPZGeoMesh;
  malha->SetName("Malha geometrica Curso Floripa");

  //Cria e insere um elemento na malha geom�rica criada acima
  UmElemento3D(*malha);

  // Os elementos da malha original podem ser refinados
  TPZVec<int> csub(0);
  TPZManVector<TPZGeoEl *> pv(4);  //vetor auxiliar que receber�os subelementos
//  int n1=1,level=0;
  cout << "\nDividir ate nivel ? ";
  int resp = 0;
//  cin >> resp;

//  int nelc = malha->ElementVec().NElements(); //nmero de elementos da malha
  int el;  //iterador

  TPZGeoEl *cpel;  //ponteiro auxiliar para um elemento geom�rico
  for(el=0;el<malha->ElementVec().NElements();el++) { //loop sobre todos os elementos da malha
    cpel = malha->ElementVec()[el];    // obt� o ponteiro para o elemento da malha de
                                       // posi�o igual ao do iterador atual
    if(cpel && cpel->Level() < resp){  // verifica se o ponteiro n� �nulo e se
                                       // se o n�el de refinamento �menor que o n�el desejado
      cpel->Divide(pv);                // Executa a divis� do elemento
    }
  }

  //solicita ao usu�io uma ordem p de interpolacao padr� para a malha
  int ord = 1;
  cout << "Entre ordem 1,2,3,4,5 : ";
//  cin >> ord;

  //Define a ordem p de cria�o de todo elemento computacional como sendo ord
//  TPZCompEl::gOrder = ord;
  TPZCompEl::SetgOrder(ord);

  // Constru�o da malha computacional. Aqui a malha ter�apenas uma refer�cia
  // para a malha geom�rica...
  TPZCompMesh *malhacomp = new TPZCompMesh(malha);

  //Pode-se dar um nome a malha para efeito de impress� de resultados...
  malhacomp->SetName("Malha Computacional :  Connects e Elementos");

  // Define-se o problema a resolver em cada parte do dom�io / definem-se os materiais
  InicializarMaterial(*malhacomp);

//  string filename("CursoAltoDes.in");
//  string filename("matriz1089.in");
//  string filename("matriz12.in");
//  string filename("matriz16000.in");
//  string filename("matriz260000.in");
  TPZFMatrix rhs,sol;
  //Define-se o espa� de aproxima�o
  malhacomp->AutoBuild();

  //Ajusta o espa� levando em conta as condi�es de contorno
  malhacomp->AdjustBoundaryElements();

  //Inicializa-se a estrutura de dados
  malhacomp->InitializeBlock();


  //Verifica�o e melhoria da banda da matriz de rigidez (Opcional)
  TPZFMatrix before, after;
  malhacomp->ComputeFillIn(100,before);

  // Cria-se um objeto an�ise. O nico par�etro �a malha computacional
  // com toda a sua estrutura de dados j�preparada. O PZ tentar�otimizar
  // a banda aqui.
  TPZAnalysis an(malhacomp);

  // c�culo da banda ap� o ajuste
  malhacomp->ComputeFillIn(100,after);

  //P� processamento da matriz de rigidez antes e depois da otimiza�o da banda
  VisualMatrix(before,"before.dx");
  VisualMatrix(after,"after.dx");

  //Escolha do padr� de armazenamento. O par�etro de entrada �a malha computacional
  //TPZSkylineStructMatrix strmat(malhacomp);
  TPZSkylineStructMatrix strmat (malhacomp);       //skyline em paralelo (multthread)
  //TPZParFrontStructMatrix<TPZFrontSym> strmat(malhacomp);

  //Define-se o padr� de armazenamento para a an�ise
  an.SetStructuralMatrix(strmat);

  //Cria�o e defini�o do solver
  TPZStepSolver step;
//  step.SetDirect(ECholesky);
  //Define-se o solver para a an�ise
  an.SetSolver(step);

//  an.Assemble();
  //TPZMatrixSolver *precond = an.BuildPreconditioner(TPZAnalysis::EBlockJacobi,false);
  TPZStepSolver jac;
  jac.SetJacobi(1,0.,0);
  jac.ShareMatrix(step);
  step.SetCG(10000,jac,1.e-10,0);
  an.SetSolver(step);
//  delete precond;
  std::ofstream out("output.txt");
  an.Assemble();
  an.Rhs().Print("Right hand side",out);

  //Processamento
  cout << "Number of equations " << malhacomp->NEquations() << endl;
  an.Run();

  an.Solution().Print("Solution obtained",out);

  // Posprocessamento
  TPZVec<std::string> scalnames(1);
  scalnames[0] = "state";    //nome das vari�eis que se quer p�-processar
  TPZVec<std::string> vecnames(0);
  std::string plotfile =  "termica.dx"; //nome do arquivo de resultados
  //Define-se para a an�ise as vari�eis a p�-processar
  an.DefineGraphMesh(2, scalnames, vecnames, plotfile);
  //Executa os c�culos para gera�o dos resultados de p�-processamento
  an.PostProcess(2);
  //   an.DefineGraphMesh(2, scalnames, vecnames, pltfile);
  //   an.PostProcess(2);

  //limpar a mem�ia alocada dinamicamente.
  delete malhacomp;
  delete malha;
  return 0;
}

void UmElemento(TPZGeoMesh &malha) {
  //Para efeito de teste ser�criado um elemento quadrilateral de
  //comprimento 1, com as seguintes coordenadas (x,y,z)
  // - canto inferior esquerdo: (0,0,0);
  // - canto inferior direito : (0,1,0);
  // - canto superior direito : (1,1,0);
  // - canto superior esquerdo: (0,1,0);
  double coordstore[4][3] = {{0.,0.,0.},{1.,0.,0.},{1.,1.,0.},{0.,1.,0.}};

  // criar os quatro n� geom�ricos
  int i,j;                              //iteradores
  TPZVec<REAL> coord(3,0.);             //vetor auxiliar para armazenar uma coordenada
  malha.NodeVec().Resize(4);

  for(i=0; i<4; i++) {                  //loop sobre o nmero de n� da malha

    // inicializar as coordenadas do no em um vetor do tipo pz
    for (j=0; j<3; j++) coord[j] = coordstore[i][j];
    // identificar um espa� no vetor da malha geom�rica
    // onde podemos armazenar o objeto n�a criar
//    int nodeindex = malha.NodeVec ().AllocateNewElement ();

    // criar um n�geom�rico e inser�lo na posi�o
    // alocada no vetor de n� da malha geom�rica
    malha.NodeVec ()[i].Initialize (i,coord,malha);

  }

  // Cria�o de um elemento geom�rico
  // inicializar os �dices dos n� do elemento
  TPZVec<int> indices(4);
  for(i=0; i<4; i++) indices[i] = i; //loop sobre o nmero de n� do elemento
  //no caso s�h�quatro n� e eles foram criados na ordem correta para o
  //elemento em quest�. A ordem dos n� deve seguir um padr� pr�estabelecido

  // O pr�rio construtor vai inserir o elemento na malha
  // os par�etros de cria�o do elemento s� os seguintes:
  // 1) Tipo geom�rico do elemento
  // - EPoint            element 0D - type point        -  associated index 0
  // - EOned             element 1D - type oned         -  associated index 1
  // - ETriangle         element 2D - type triangle     -  associated index 2
  // - EQuadrilateral    element 2D - type quad         -  associated index 3
  // - ETetraedro        element 3D - type tetraedro    -  associated index 4
  // - EPiramide         element 3D - type piramide     -  associated index 5
  // - EPrisma           element 3D - type prisma       -  associated index 6
  // - ECube             element 3D - type cube         -  associated index 7
  // 2) Vetor de conectividades dos elementos (o nmero de n� deve ser
  //    compat�el com o nmero de n� do tipo do elemento
  // 3) �dice do material que ser�associado ao elemento
  // 4) Vari�el onde ser�retornado o �dice do elemento criado no vetor
  //    de elementos da malha geom�rica
  // TPZGeoEl *gel = new TPZGeoElQ2d(0,indices,1,malha);  //forma antiga.
  int index;
  TPZGeoEl *gel = malha.CreateGeoElement(EQuadrilateral,indices,1,index);

  //Gerar as estruturas de dados de conectividade e vizinhan�
  malha.BuildConnectivity ();

  // Identificar onde ser� inseridas condi�es de contorno.
  // Uma condi�o de contorno �aplicada a um lado (parte do contorno)
  // de um elemento. Este objeto ir�inserir-se automaticamente na malha
  // Os par�etros s� os seguintes:
  // 1) elemento onde ser�aplicada a condi�o de contorno
  // 2) lado do elemento onde ser�inserida a condi�o de contorno
  // 3) identificador da condi�o de contorno
  // 4) refer�cia para a malha geom�rica.
  TPZGeoElBC(gel,4,-4);
}

void UmElemento3D(TPZGeoMesh &malha) {
  //Para efeito de teste ser�criado um elemento quadrilateral de
  //comprimento 1, com as seguintes coordenadas (x,y,z)
  // - canto inferior esquerdo: (0,0,0);
  // - canto inferior direito : (0,1,0);
  // - canto superior direito : (1,1,0);
  // - canto superior esquerdo: (0,1,0);
  double coordstore[8][3] = {{0.,0.,0.},{1.,0.,0.},{1.,1.,0.},{0.,1.,0.},
	{0.,0.,1.},{1.,0.,1.},{1.,1.,1.},{0.,1.,1.}};

  // criar os quatro n� geom�ricos
  int i,j;                              //iteradores
  TPZVec<REAL> coord(3,0.);             //vetor auxiliar para armazenar uma coordenada
  malha.NodeVec().Resize(8);

  for(i=0; i<8; i++) {                  //loop sobre o nmero de n� da malha

    // inicializar as coordenadas do no em um vetor do tipo pz
    for (j=0; j<3; j++) coord[j] = coordstore[i][j];
    // identificar um espa� no vetor da malha geom�rica
    // onde podemos armazenar o objeto n�a criar
//    int nodeindex = malha.NodeVec ().AllocateNewElement ();

    // criar um n�geom�rico e inser�lo na posi�o
    // alocada no vetor de n� da malha geom�rica
    malha.NodeVec ()[i].Initialize (i,coord,malha);

  }

  // Cria�o de um elemento geom�rico
  // inicializar os �dices dos n� do elemento
  TPZVec<int> indices(8);
  for(i=0; i<8; i++) indices[i] = i; //loop sobre o nmero de n� do elemento
  //no caso s�h�quatro n� e eles foram criados na ordem correta para o
  //elemento em quest�. A ordem dos n� deve seguir um padr� pr�estabelecido

  // O pr�rio construtor vai inserir o elemento na malha
  // os par�etros de cria�o do elemento s� os seguintes:
  // 1) Tipo geom�rico do elemento
  // - EPoint            element 0D - type point        -  associated index 0
  // - EOned             element 1D - type oned         -  associated index 1
  // - ETriangle         element 2D - type triangle     -  associated index 2
  // - EQuadrilateral    element 2D - type quad         -  associated index 3
  // - ETetraedro        element 3D - type tetraedro    -  associated index 4
  // - EPiramide         element 3D - type piramide     -  associated index 5
  // - EPrisma           element 3D - type prisma       -  associated index 6
  // - ECube             element 3D - type cube         -  associated index 7
  // 2) Vetor de conectividades dos elementos (o nmero de n� deve ser
  //    compat�el com o nmero de n� do tipo do elemento
  // 3) �dice do material que ser�associado ao elemento
  // 4) Vari�el onde ser�retornado o �dice do elemento criado no vetor
  //    de elementos da malha geom�rica
  // TPZGeoEl *gel = new TPZGeoElQ2d(0,indices,1,malha);  //forma antiga.
  int index;
  TPZGeoEl *gel = malha.CreateGeoElement(ECube,indices,1,index);

  //Gerar as estruturas de dados de conectividade e vizinhan�
  malha.BuildConnectivity ();

  // Identificar onde ser� inseridas condi�es de contorno.
  // Uma condi�o de contorno �aplicada a um lado (parte do contorno)
  // de um elemento. Este objeto ir�inserir-se automaticamente na malha
  // Os par�etros s� os seguintes:
  // 1) elemento onde ser�aplicada a condi�o de contorno
  // 2) lado do elemento onde ser�inserida a condi�o de contorno
  // 3) identificador da condi�o de contorno
  // 4) refer�cia para a malha geom�rica.
  TPZGeoElBC(gel,0,-1);
  TPZGeoElBC(gel,1,-2);
  TPZGeoElBC(gel,2,-3);
  TPZGeoElBC(gel,24,-4);
  TPZGeoElBC(gel,22,-5);

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
// os par�etros utilizados s� o nome do arquivo e uma
// refer�cia para a malha geom�rica que receber�os dados
void LerMalha(char *nome, TPZGeoMesh &grid) {
  //Cria�o de um objeto arquivo de entrada de dados;
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
		TPZGeoElBC(el,side,-mat);
	}
	grid.BuildConnectivity();

	return;
}

void InicializarMaterial(TPZCompMesh &cmesh) {
  //O par�etro b�ico para a cria�o de um material �o
  // identificador. Esse deve corresponder aos �dices utilizados
  // na cria�o dos elementos geom�ricos.
  // No restante, o material define a equa�o diferencial a ser
  // resolvida...

  int dim = cmesh.Reference()->ElementVec()[0]->Dimension();
  TPZMat2dLin *meumatp = new TPZMat2dLin(2);
  TPZMatPoisson3d *poismatp = new TPZMatPoisson3d(3,dim);
  TPZVec<REAL> force(3,0.);
  TPZElasticity3D *elasmatp = new TPZElasticity3D(1,1000.,0.,force);
  TPZVec<REAL> convdir(3,0.);
  poismatp->SetParameters(1.,0.,convdir);
  poismatp->SetInternalFlux(1.);
  //poismat->SetForcingFunction(forcingfunction);

  //Cada material tem par�etros de inicializa�o pr�rios, assim
  //deve-se consultar a documenta�o para verificar como definir os
  //par�etros. No caso em quest� o material requer tr� matrizes
  //e uma fun�o de c�culo tamb� �fornecida
  TPZFMatrix xk(1,1,1.),xc(1,2,0.),xf(1,1,1.);
  meumatp->SetMaterial (xk,xc,xf);
//  meumat->SetForcingFunction(forcingfunction);

  //Ap� a cria�o do material este dever ser inserido na estrutura
  //de dados da malha computacional
  TPZAutoPointer<TPZMaterial> meumat(meumatp);
  TPZAutoPointer<TPZMaterial> poismat(poismatp);
  TPZAutoPointer<TPZMaterial> elasmat(elasmatp);
  cmesh.InsertMaterialObject(meumat);
  cmesh.InsertMaterialObject(poismat);
  cmesh.InsertMaterialObject(elasmat);

//  TPZMaterial *atual = meumat;
//  TPZMaterial *atual = poismat;
  TPZAutoPointer<TPZMaterial> atual = elasmat;

  // inserir as condi�es de contorno
  // Uma condi�o de contorno pode ser dada por duas matrizes
  // relacionadas a rigidez e ao vetor de carga. As dimens�s dessas
  // matrizes dependem da dimens� e do nmero de vari�eis de estado
  // do problema. O tipo de condi�o de contorno �definido por um
  // identificador, onde:
  //  - 0 = Dirichlet
  //  - 1 = Neumann
  //  - 2 = Mista
  int nstate = atual->NStateVariables();
  TPZFMatrix val1(nstate,nstate,0.),val2(nstate,1,0.);

  //Os par�etros necess�ios �cria�o de uma condi�o de contorno s�:
  // 1) Identificador da condi�o de contorno (lembre-se que uma BC �como um material)
  // 2) O Tipo da BC : Dirichlet, Neumann ou Mista
  // 3) Os valores da BC

  if(nstate == 3)
  {
    val2(1,0) = -1.;
  }
  TPZAutoPointer<TPZMaterial> bnd = atual->CreateBC (atual,-4,1,val1,val2);
  cmesh.InsertMaterialObject(bnd);
  if(nstate == 3)
  {
    val2(1,0) = 1.;
  }
  bnd = atual->CreateBC (atual,-5,1,val1,val2);
  val2.Zero();

  //Da mesma forma que para os materiais, ap� sua cria�o �necess�io
  //a sua inser�o na estrutura de dados da malha computacional
  cmesh.InsertMaterialObject(bnd);

  //cria�o e inser�o de outras BC's
  bnd = atual->CreateBC (atual,-1,0,val1,val2);
  cmesh.InsertMaterialObject(bnd);

  if(nstate == 3)
  {
    val1(1,1) = 1000.;
    val1(2,2) = 1000.;
  }
  bnd = atual->CreateBC (atual,-3,2,val1,val2);
  cmesh.InsertMaterialObject(bnd);
  val1.Zero();
  val2.Zero();
  if(nstate == 3)
  {
    val1(2,2) = 1000.;
  }
  bnd = atual->CreateBC (atual,-2,2,val1,val2);
  cmesh.InsertMaterialObject(bnd);
}

