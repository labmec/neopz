///#include "pzmetis.h"
//#include "pztrnsform.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzfmatrix.h"
#include "pzelgc3d.h"
#include "pzbndcond.h"
//#include "pztempmat.h"
#include "pzcompel.h"
#include "pzanalysis.h"
#include "pzskylstrmatrix.h"
#include "pzskylmat.h"
#include "pzstepsolver.h"
#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzstack.h"
#include "pzvec.h"
#include "pzsolve.h"
#include "pzelgpoint.h"
#include "pzelg1d.h"
#include "pzelgq2d.h"
#include "pzelgt2d.h"
#include "pzelct2d.h"
#include "pzelcc3d.h"
#include "pzelgt3d.h"
#include "pzelct3d.h"
#include "pzelgpi3d.h"
#include "pzelcpi3d.h"
#include "pzelgpr3d.h"
#include "pzelcpr3d.h"
#include "pzelmat.h"
#include "pzelasmat.h"
#include "pzmattest.h"
#include "pzmat1dlin.h"
#include "pzmat2dlin.h"
#include "pzpoisson3d.h"
#include "TPZConservationLaw.h"
#include "TPZDiffusionConsLaw.h"
#include "TPZConsLawTest.h"
//#include "TPZRefPattern.h"
#include "TPZCompElDisc.h"
#include "TPZShapeDisc.h"
#include "TPZInterfaceEl.h"
//#include "TPZJacobMat.h"
//#include "TPZJacobStrMatrix.h"
#include "pzdxmesh.h"


#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <ostream>
//#define NOTDEBUG

/* static double hexaedro[8][3] = { {0.,0.,0.},{1.,0.,0.},{1.,1.,0.},{0.,1.,0.}, */
/* 		             {0.,0.,1.},{1.,0.,1.},{1.,1.,1.},{0.,1.,1.} };//0 1 2 3 4 5 6 7 */

//static double triangulo[3][3] = { {0.,0.,-1.},{2.,0.,0.},{0.,0.,1.} };//0 1 2
static double triangulo[3][3] = { {0.,0.,0.},{1.,0.,0.},{0.,1.,0.} };//0 1 2

static double quadrilatero[4][3] = { {0.,0.,0.},{1.,0.,0.},{1.,1.,0.},{0.,1.,0.} };//0 1 2 3

static double tetraedro[4][3] = { {0.,0.,0.},{1.,0.,0.},{0.,1.,0.},{0.,0.,1.} };//0 1 2 3 4

// static double prisma[6][3] = { {0.,0.,0.},{1.,0.,0.},{0.,1.,0.},
// 		             {0.,0.,1.},{1.,0.,1.},{0.,1.,1.} };//0 1 2 3 4 5

static double prisma[6][3] = { {0.,0.,-1.},{1.,0.,-1.},{0.,1.,-1.},
			       {0.,0.,1.},{1.,0.,1.},{0.,1.,1.} };//0 1 2 3 4 5

static double piramide[5][3] = { {0.,0.,0.},{1.,0.,0.},{1.,1.,0.},{0.,1.,0.},{0.5,0.5,1.} };//0 1 2 3 4

static double hexaedro[8][3] = { {0.,0.,0.},{1.,0.,0.},{1.,1.,0.},{0.,1.,0.},
				 {0.,0.,1.},{1.,0.,1.},{1.,1.,1.},{0.,1.,1.} };//0 1 2 3 4 5 6 7

static double linha[2][3] = { {0.,0.,0.},{1.,0.,0.}};

void CriacaoDeNos(int nnodes,double lista[20][3]);
TPZMaterial *Linha(int ordem);
TPZMaterial *Triangulo(int ordem);
TPZMaterial *Quadrilatero(int ordem);
TPZMaterial *QuadrilateroNovo(int ordem);
TPZMaterial *Hexaedro(int ordem);
TPZMaterial *HexaedroNovo(int ordem);
void Tetraedro();
void Piramide();
void Prisma();
void ContagemDeElementos();
void CycleRefinements(TPZCompMesh &cm, int numcycles, int minel, int maxel, ofstream &out);
int IsGroup(int ind, TPZCompEl *cel, TPZVec<int> &indgroup);
void FileNB(TPZGeoMesh &gmesh,ostream &out);
void BCZero(TPZVec<REAL> &x,TPZVec<REAL> &result);
void CC_Poisson(TPZMaterial *mater);
void CC_Projecao(TPZMaterial *mater);
void Function(TPZVec<REAL> &x,TPZVec<REAL> &result);
void BCPoisson(TPZVec<REAL> &x,TPZVec<REAL> &result);
void Solution(TPZVec<REAL> &x,TPZVec<REAL> &result,TPZFMatrix &deriv);
void Divisao (TPZCompMesh *cmesh,int key);
void PostProcess(TPZGeoMesh &gmesh,ostream &out);
//static TPZRefPattern hexaedro9("hexaedro9subs.in");
//static TPZRefPattern hexaedro8("hexaedro8.in");//hexaedro mestre 8 sub-elementos
//static TPZRefPattern hexaedro2("hexaedro2.in");
//static TPZRefPattern quad4("quadrilatero4.in");//quadrilatero mestre 4 sub-elementos*/
//static TPZRefPattern linha2("linha2.in");//aresta mestre 2 sub-elementos*/
static TPZGeoMesh *gmesh = new TPZGeoMesh;
static TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
void TestShapesDescontinous();
//static int problemtype;
static TPZVec<REAL> x0(3,0.);
void NivelDivide(TPZCompMesh *cmesh);
static int grau = 0;
static int nivel = 0;
static int problem=0;
static REAL r=0.2,pi = 2.0*asin(1.0);
static clock_t start,end;//,begin,ttot=0;
void CoutTime(clock_t &start);
TPZMaterial *mat=0;
REAL CFL = -1;


int main() {


  cout << "\nProblema (0/1/9/10): ";
  cin >> problem;
  x0[0] = 0.2;//problem 10 em 1d
  if(problem == 9){
    x0[0] = 0.2;
    x0[1] = 0.2;
    x0[2] = 0.5;
  }
  if(problem == 10){
    cout << "\nEntre 3 coordenadas do ponto centro do chapeu\n";
    cin >> x0[0] >> x0[1] >> x0[2];
  }
  int tipo;
  cout << "\nElemento [0:linha][1:triangulo][2:quadrilatero][3:tetraedro][4:piramide][5:prisma][6:hexaedro]\n\t";
  cin >> tipo;
  cout << "\nGrau do espaco de interpolacao -> 0,1,2,3,...";
  cin >> grau;
  TPZCompElDisc::gDegree = grau;

  if(tipo==0) mat = Linha(grau);

  if(tipo==1) mat = Triangulo(grau);

  if(tipo==2) mat = QuadrilateroNovo(grau);

  if(tipo==6) mat = HexaedroNovo(grau);

  if(tipo==3) Tetraedro();

  if(tipo==4) Piramide();

  if(tipo==5) Prisma();

  ofstream outgm("mesh.out");

  if(0){//cyclerefinements
    int minel = 10,maxel = 20,numcyc = 2;//500,1000
    cout << "\nEntre minel : ";
    cin >> minel;
    cout << "\nEntre maxel : ";
    cin >> maxel;
    cout << "\nEntre numcycles : ";
    cin >> numcyc;
    CycleRefinements(*cmesh,numcyc,minel,maxel,outgm);//maxel=150, minel = 2, numcyc=50
  }

  if(0){
    cout << "\nmain::Divisao manual\n";
    Divisao (cmesh,1);
  }

  if(1){
    cout << "\ndescontinuo.c::main verificando a consistencia da malha de interfaces\t";
    if(TPZInterfaceElement::main(*cmesh)){
      cout << "->\tOK!";
    } else {
      cout << "->\tPROBLEMAS COM INTERFACES\n\n";
      //return 0;
    }
    //ContagemDeElementos();
  }

  if(0){
    cout << "\nmain::Imprime malhas\n";
    gmesh->Print(outgm);
    cmesh->Print(outgm);
    outgm.flush();
  }

  if(1){
    start = clock();
    NivelDivide(cmesh);
    CoutTime(start);
    if(0){
      gmesh->Print(outgm);
      cmesh->Print(outgm);
      outgm.flush();
    }
  }

  if(1){
    cout << "\nmain::Ajuste no contorno\n";
    cmesh->AdjustBoundaryElements();
    if(0){
      cout << "\nmain::Imprime malhas\n";
      gmesh->Print(outgm);
      cmesh->Print(outgm);
      outgm.flush();
    }
  }

  if(1){
    TPZAnalysis an(cmesh,outgm);
    if(1){//Analysis
      cout << "\nmain::Resolve o sistema\n";
      //int nel2d,nel3d,cap1,cap2;
      //TPZFMatrix bloco = cmesh->BlockJacob(cap1,cap2,nel2d,nel3d);
      //TPZJacobStructMatrix stiff(cap1,cap2,nel2d,nel3d,bloco);
      //TPZJacobStructMatrix stiff(cmesh);
      TPZSkylineStructMatrix stiff(cmesh);
      an.SetStructuralMatrix(stiff);
      an.Solution().Zero();
      TPZStepSolver solver;//,precond;
      solver.SetDirect(ELDLt);//ELU, ECholesky
//       int numiterations=100,numvectors=10,FromCurrent=1000;
//       REAL tol = 1.e-8;
//       if(1){
// 	cout << "main:: Entre num iter, num vetor,num max iter, tol : \n";
// 	cin >> numiterations >> numvectors >> FromCurrent >> tol;
//       }
//       precond.SetDirect(ELU);
//       solver.SetGMRES(numiterations,numvectors,precond,tol,FromCurrent);
      an.SetSolver(solver);
      if(1){
	int numiter=20;
	REAL tol;
	cout << "\nNumero de iteracoes? : ";
	cin >> numiter;
	tol = 1.0e15;// = norma da solu��o inicial + epsilon
	cout << "\nTolerancia? : " << tol << "\n";
	//cin >> tol;
	an.SetExact(Solution);
	int marcha=0,resolution=0;
	cout << "main:: Parametro marcha : \n";
	cin >> marcha;
	//marcha = 1;
	cout << "main:: Parametro resolution : \n";
	//cin >> resolution;
	resolution = 0;
	an.IterativeProcessTest(outgm,tol,numiter,mat,marcha,resolution);
	if(0) PostProcess(*gmesh,outgm);
	if(mat->Dimension() == 1){
	  ofstream plot("APLOT.nb");
	  FileNB(*gmesh,plot);
	}
      }
      if(0)
	{//pos-processamento
	  TPZVec<char *> scalar(1),vector(0);
	  scalar[0] = "Solution";
	  //vector[0] = "Solution";//"Tension9";
	  char plot[] = "conslaw.dx";
	  an.DefineGraphMesh(2,scalar,vector,plot);
	  an.PostProcess(0);
	}
    }

    //ofstream outan("analysis.out");
    if(0) an.Print("FEM SOLUTION ",outgm);

    REAL estimate=0.,H1_error=0.,L2_error=0.;

    if(0){//erro da aproximacao FEM
      TPZManVector<REAL> flux(9);
      cmesh->EvaluateError(Solution,H1_error,L2_error,estimate);
      outgm << "\n\n\n";
      outgm << "L2 Error  : " << L2_error << endl;
      outgm << "Semi Norm : " << estimate << endl;
      outgm << "Energy    : " << H1_error << endl;
      outgm.flush();
      cout << "L2 Error  : " << L2_error << endl;
      cout << "Semi Norm : " << estimate << endl;
      cout << "Energy    : " << H1_error << endl;
    }
    if(0){//pos-processamento local
      PostProcess(*gmesh,outgm);
      outgm.flush();
    }

    ContagemDeElementos();
  }//if(0/1)

  outgm.close();
  if(cmesh) delete cmesh;
  if(gmesh) delete gmesh;
  //AvisoAudioVisual();
  return 0;
}

void CriacaoDeNos(int nnodes,double lista[20][3]){

   gmesh->NodeVec().Resize(nnodes);
   TPZVec<REAL> coord(3);
   int i;
   for(i=0;i<nnodes;i++){
     coord[0] = lista[i][0];
     coord[1] = lista[i][1];
     coord[2] = lista[i][2];
     gmesh->NodeVec()[i].Initialize(coord,*gmesh);
  }
}

TPZMaterial *Linha(int ordem){

  CriacaoDeNos(2,linha);
  //elemento de volume
  TPZVec<int> nodes;
  nodes.Resize(2);
  nodes[0] = 0;
  nodes[1] = 1;
  TPZGeoEl1d *elg1d = new TPZGeoEl1d(nodes,1,*gmesh);
  TPZGeoEl1d::SetCreateFunction(TPZCompElDisc::Create1dDisc);//de volume
  TPZGeoElPoint::SetCreateFunction(TPZCompElDisc::CreatePointDisc);//n�o de interface
  int interfdim = 0;
  TPZCompElDisc::gInterfaceDimension = interfdim;
  gmesh->BuildConnectivity();
  TPZVec<REAL> B(1);
  B[0] = 1.0;
  int nummat = 1;
  int dif_artif = 1;//SUPG
  cout << "\nmain::Divisao Nivel final da malha ? : ";
  cin >> nivel;
  CFL = ( 1./(2.0*(REAL)grau+1.0) );
  REAL delta_x = 1.0 / pow(2.0,(REAL)nivel);//tamanho do elemento
  REAL delta = ( (10./3.)*CFL*CFL - (2./3.)*CFL + 1./10. );
  REAL maxflinha = 1.0;
  REAL delta_t = CFL*delta_x / maxflinha;//delta_t � <= que este valor
  cout << "\nMaximo de f' = " << maxflinha
       << "\nCFL = 1/(2p+1) = " << CFL
       << "\nDominio [0,1]"
       << "\nDelta x = " << delta_x
       << "\ndelta t = " << delta_t
       << "\ndelta otimo = " << delta << endl;
  int dim = 1;
  //int integ = 0;
  //cout << "\nGrau de integracao do problema : ";
  //cin >> integ;
  int nao;
  cout << "main::Linha modificar CFL? ";
  cin >> nao;
  if(nao){
    cout << "main::SetDeltaTime novo CFL -> ";
    cin >> CFL;
    TPZDiffusionConsLaw::fCFL = CFL;
    delta_t = CFL*delta_x/maxflinha;
    cout << "main::Linha modificar delta? ";
    cin >> nao;
    if(nao){
      delta = ( (10./3.)*CFL*CFL - (2./3.)*CFL + 1./10. );
    }
    cout << "\nNovo delta t = " << delta_t
	 << "\nNovo delta   = " << delta << endl;
  }
  TPZMaterial *mater = new TPZConsLawTest(nummat,B,dif_artif,delta_t,dim,delta);
  dynamic_cast<TPZConservationLaw *>(mater)->SetTimeStep(delta_t);
  //((TPZConservationLaw *)mat)->SetIntegDegree(grau);
  mater->SetForcingFunction(Function);
  cmesh->InsertMaterialObject(mater);
  //condi��es de contorno
  TPZBndCond *bc;
  TPZFMatrix val1(1,1,0.),val2(1,1,0.);

  //CC no canto -1
  val1.Zero();
  val2.Zero();
  if(problem==0) val2(0,0) = 1.;
  TPZGeoElBC(elg1d,0,-1,*gmesh);
  bc = mater->CreateBC(-1,3,val1,val2);
  cmesh->InsertMaterialObject(bc);

  //CC no canto 1
  TPZGeoElBC(elg1d,1,-2,*gmesh);
  bc = mater->CreateBC(-2,3,val1,val2);
  //bc->SetForcingFunction(BCZero);
  cmesh->InsertMaterialObject(bc);

  cout << endl;
  cmesh->AutoBuild();

  return mater;
}

void Function(TPZVec<REAL> &x,TPZVec<REAL> &result){

  if(problem==0) result[0] = 1.0;

  if(problem==1) result[0] = x[0];

  if(problem==9 || problem==10){
    int i;//,nel = x.NElements();
    REAL norma = 0.;
    for(i=0;i<3;i++) norma += (x[i]-x0[i])*(x[i]-x0[i]);
    norma = sqrt(norma);//norma(x-x0)
    if(norma <= r){
      REAL teta = norma*pi/r;
      result[0] = 1.0+cos(teta);
      return;
    }
    result[0] = 0.0;
  }
}

void BCZero(TPZVec<REAL> &x,TPZVec<REAL> &result){

  int nel = result.NElements(),i;
  for(i=0;i<nel;i++) result[i] = 0.0;
}

TPZMaterial *Quadrilatero(int ordem){

  CriacaoDeNos(4,quadrilatero);
  //elemento de volume
  TPZVec<int> nodes;
  nodes.Resize(4);
  nodes[0] = 0;
  nodes[1] = 1;
  nodes[2] = 2;
  nodes[3] = 3;
  TPZGeoElQ2d *elgq2d = new TPZGeoElQ2d(nodes,1,*gmesh);
  TPZGeoElQ2d::SetCreateFunction(TPZCompElDisc::CreateQ2Disc);//de volume
  TPZGeoEl1d::SetCreateFunction(TPZCompElDisc::Create1dDisc);//n�o de interface
  int interfdim = 1;
  TPZCompElDisc::gInterfaceDimension = interfdim;
  gmesh->BuildConnectivity();
  TPZVec<REAL> B(2);
  int test = 0;
  if(!problem){
    B[0] = 1.0;//1d e 2d
    B[1] = 0.0;//test=0
  }
  if(problem==1){
    B[0] = 1.0;//1d e 2d
    B[1] = 0.0;//test=0
  }
  if(problem==9){
    B[0] = 1.0;//1d e 2d
    B[1] = 0.0;//test=0
  }
  if(problem==10){
    test = 1;
  }
  int nummat = 1;
  int dif_artif = 1;//SUPG
  cout << "\nmain::Divisao Nivel final da malha ? : ";
  cin >> nivel;
  CFL = ( 1./(2.0*(REAL)grau+1.0) );
  REAL delta_x = ( 1.0 / pow(2.0,(REAL)nivel) );//tamanho do elemento
  REAL delta = ( (10./3.)*CFL*CFL - (2./3.)*CFL + 1./10. );
  REAL maxflinha = 1.0;
  REAL delta_t = CFL*delta_x / maxflinha;//delta_t � <= que este valor
  cout << "\nMaximo de f' = " << maxflinha
       << "\nCFL = " << CFL
       << "\nDominio [0,1]x[0,1]"
       << "\nDelta x = " << delta_x
       << "\ndelta t = " << delta_t
       << "\ndelta = " << delta << endl;
  int dim = 2;
  int nao;
  cout << "main::Linha modificar CFL? ";
  cin >> nao;
  if(nao){
    cout << "main::SetDeltaTime novo CFL -> ";
    cin >> CFL;
    TPZDiffusionConsLaw::fCFL = CFL;
    delta_t = CFL*delta_x/maxflinha;
    cout << "main::Linha modificar delta? ";
    cin >> nao;
    if(nao){
      delta = ( (10./3.)*CFL*CFL - (2./3.)*CFL + 1./10. );
    }
    cout << "\nNovo delta t = " << delta_t
	 << "\nNovo delta   = " << delta << endl;
  }
  TPZMaterial *mater = new TPZConsLawTest(nummat,B,dif_artif,delta_t,dim,delta);
  dynamic_cast<TPZConservationLaw *>(mater)->SetTimeStep(delta_t);
  //((TPZConservationLaw *)mater)->SetIntegDegree(grau);
  mater->SetForcingFunction(Function);
  cmesh->InsertMaterialObject(mater);
  //condi��es de contorno
  TPZBndCond *bc;
  //  REAL big = 1.e12;
  TPZFMatrix val1(1,1,0.),val2(1,1,0.);

  //a �nica fun��o atual destes elemento
  //� a cria��o dos elementos de interface
  //CC na aresta 4
  val1.Zero();//para uso CC tipo 2
  val2.Zero();

  if(problem==0) val2(0,0) = 1.;
  TPZGeoElBC(elgq2d,4,-1,*gmesh);
  bc = mater->CreateBC(-1,1,val1,val2);
  //bc->SetForcingFunction(Function);
  cmesh->InsertMaterialObject(bc);

  //CC na aresta 5
  TPZGeoElBC(elgq2d,5,-2,*gmesh);
  bc = mater->CreateBC(-2,1,val1,val2);
  //bc->SetForcingFunction(Function);
  cmesh->InsertMaterialObject(bc);

  //CC na aresta 6
  TPZGeoElBC(elgq2d,6,-3,*gmesh);
  bc = mater->CreateBC(-3,1,val1,val2);
  //bc->SetForcingFunction(Function);
  cmesh->InsertMaterialObject(bc);

  //CC na aresta 7
  TPZGeoElBC(elgq2d,7,-4,*gmesh);
  bc = mater->CreateBC(-4,1,val1,val2);
  //bc->SetForcingFunction(Function);
  cmesh->InsertMaterialObject(bc);

  cout << endl;
  cmesh->AutoBuild();

  return mater;
}

TPZMaterial *Hexaedro(int ordem){

  CriacaoDeNos(8,hexaedro);
  //elemento de volume
  TPZVec<int> nodes;
  nodes.Resize(8);
  nodes[0] = 0;
  nodes[1] = 1;
  nodes[2] = 2;
  nodes[3] = 3;
  nodes[4] = 4;
  nodes[5] = 5;
  nodes[6] = 6;
  nodes[7] = 7;
  TPZGeoElC3d *elgc3d = new TPZGeoElC3d(nodes,1,*gmesh);
  //construtor descont�nuo
  TPZGeoElC3d::SetCreateFunction(TPZCompElDisc::CreateC3Disc);//de volume
  TPZGeoElQ2d::SetCreateFunction(TPZCompElDisc::CreateQ2Disc);//n�o de interface
  int interfdim = 2;
  TPZCompElDisc::gInterfaceDimension = interfdim;
  gmesh->BuildConnectivity();
  TPZVec<REAL> B(3);
  int test = 0;
  if(!problem){
    B[0] = 1.0;
    B[1] = 0.0;
    B[2] = 0.0;
  }
  if(problem==1){
    B[0] = 1.0;
    B[1] = 0.0;
    B[2] = 0.0;
  }
  if(problem==9){
    B[0] = 1.0;
    B[1] = 0.0;
    B[2] = 0.0;
  }
  if(problem==10){
    test = 2;
  }
  int nummat = 1;
  int dif_artif = 1;//SUPG
  cout << "\nmain::Divisao Nivel final da malha ? : ";
  cin >> nivel;
  REAL cfl = .8;//( 1./(2.0*(REAL)grau+1.0) );
  REAL delta_x = ( 1.0 / pow(2.0,(REAL)nivel) );//tamanho do elemento
  REAL delta = 2.5*delta_x;//( (10./3.)*cfl*cfl - (2./3.)*cfl + 1./10. );
  REAL maxflinha = 1.0;
  REAL delta_t = cfl*delta_x / maxflinha;//delta_t � <= que este valor
  //calculando novos valores
  delta_t = delta_x*cfl;
  //delta = delta_x / 2.0;
  cout << "\nMaximo de f' = " << maxflinha
       << "\nCFL = " << cfl
       << "\nDominio [0,1]x[0,1]x[0,1]"
       << "\nDelta x = " << delta_x
       << "\ndelta t = " << delta_t
       << "\ndelta = " << delta << endl;
//   cout << "main::Hexaedro Entre com o passo : ";
//   cin >> delta_t;
  int dim = 3;
//   int integ = 0;
//   cout << "\nmain::Hexaedro Grau de integracao do problema : ";
//   cin >> integ;
  TPZMaterial *mater = new TPZConsLawTest(nummat,B,dif_artif,delta_t,dim,delta,test);
  //((TPZConservationLaw *)mater)->SetIntegDegree(grau);
  mater->SetForcingFunction(Function);
  cmesh->InsertMaterialObject(mater);
  //condi��es de contorno
  TPZBndCond *bc;
  //  REAL big = 1.e12;
  TPZFMatrix val1(1,1,0.),val2(1,1,0.);

  //a �nica fun��o atual destes elemento
  //� a cria��o dos elementos de interface
  //CC na aresta 4
  val1.Zero();//para uso CC tipo 2
  val2.Zero();

  if(problem==0) val2(0,0) = 1.;
  //CC na aresta 20
  TPZGeoElBC(elgc3d,20,-1,*gmesh);
  bc = mater->CreateBC(-1,1,val1,val2);
  //bc->SetForcingFunction(Function);
  cmesh->InsertMaterialObject(bc);

  //CC na aresta 21
  TPZGeoElBC(elgc3d,21,-2,*gmesh);
  bc = mater->CreateBC(-2,1,val1,val2);
  //bc->SetForcingFunction(Function);
  cmesh->InsertMaterialObject(bc);

  //CC na aresta 22
  TPZGeoElBC(elgc3d,22,-3,*gmesh);
  bc = mater->CreateBC(-3,1,val1,val2);
  //bc->SetForcingFunction(Function);
  cmesh->InsertMaterialObject(bc);

  //CC na aresta 23
  TPZGeoElBC(elgc3d,23,-4,*gmesh);
  bc = mater->CreateBC(-4,1,val1,val2);
  //bc->SetForcingFunction(Function);
  cmesh->InsertMaterialObject(bc);

  //CC na aresta 24
  TPZGeoElBC(elgc3d,24,-5,*gmesh);
  bc = mater->CreateBC(-5,1,val1,val2);
  //bc->SetForcingFunction(Function);
  cmesh->InsertMaterialObject(bc);

  //CC na aresta 25
  TPZGeoElBC(elgc3d,25,-6,*gmesh);
  bc = mater->CreateBC(-6,1,val1,val2);
  //bc->SetForcingFunction(Function);
  cmesh->InsertMaterialObject(bc);

//   int ordem=2;
//   cout << "\nEntre ordem -> ";
//   //cin >> ordem;
//   TPZCompEl::gOrder = ordem;
  //cmesh->SetDefaultOrder(ordem);
  cout << endl;
  cmesh->AutoBuild();

  return mater;
}

TPZMaterial *Triangulo(int ordem){

  CriacaoDeNos(3,triangulo);
  //elemento de volume
  TPZVec<int> nodes;
  nodes.Resize(3);
  nodes[0] = 0;
  nodes[1] = 1;
  nodes[2] = 2;
  TPZGeoElT2d *elgt2d = new TPZGeoElT2d(nodes,1,*gmesh);
  TPZGeoElT2d::SetCreateFunction(TPZCompElDisc::CreateT2Disc);//de volume
  TPZGeoEl1d::SetCreateFunction(TPZCompElDisc::Create1dDisc);//n�o de interface
  int interfdim = 1;
  TPZCompElDisc::gInterfaceDimension = interfdim;
  gmesh->BuildConnectivity();
  TPZVec<REAL> B(2);
  B[0] = 1.0;
  B[1] = 0.0;
  int nummat = 1;
  int dif_artif = 1;//SUPG
  cout << "main::Triangulo Entre com o passo : ";
  REAL delta_t = 0.1;
  cin >> delta_t;
  int dim = 2,integ = 0;
  cout << "\nmain::Triangulo Grau de integracao do problema : " << integ << endl;
  TPZMaterial *mater = new TPZConsLawTest(nummat,B,dif_artif,delta_t,dim,integ);
  mater->SetForcingFunction(Function);
  cmesh->InsertMaterialObject(mater);
  //condi��es de contorno
  TPZBndCond *bc;
  //  REAL big = 1.e12;
  TPZFMatrix val1(1,1,0.),val2(1,1,0.);

  //a �nica �nica fun��o atual destes elemento
  //� a cria��o dos elementos de interface
  //CC na aresta 3
  val1.Zero();//para uso CC tipo 2
  val2.Zero();
  TPZGeoElBC(elgt2d,4,-1,*gmesh);
  bc = mater->CreateBC(-1,1,val1,val2);
  bc->SetForcingFunction(BCZero);
  cmesh->InsertMaterialObject(bc);

  //CC na aresta 4
  TPZGeoElBC(elgt2d,5,-2,*gmesh);
  bc = mater->CreateBC(-2,1,val1,val2);
  bc->SetForcingFunction(BCZero);
  cmesh->InsertMaterialObject(bc);

  //CC na aresta 5
  TPZGeoElBC(elgt2d,6,-3,*gmesh);
  bc = mater->CreateBC(-3,1,val1,val2);
  bc->SetForcingFunction(BCZero);
  cmesh->InsertMaterialObject(bc);

  cout << endl;
  cmesh->AutoBuild();

  return mater;
}

void Piramide(){

  CriacaoDeNos(5,piramide);
  //elemento de volume
  TPZVec<int> nodes;
  nodes.Resize(5);
  nodes[0] = 0;
  nodes[1] = 1;
  nodes[2] = 2;
  nodes[3] = 3;
  nodes[4] = 4;
  new TPZGeoElPi3d(nodes,100,*gmesh);//TPZGeoElPi3d *elg3d =
  int interfdim = 2;
  TPZCompElDisc::gInterfaceDimension = interfdim;
  TPZGeoElPi3d::SetCreateFunction(TPZCompElDisc::CreatePi3Disc);
  TPZGeoElT3d::SetCreateFunction(TPZCompElDisc::CreateT3Disc);
  TPZGeoElT2d::SetCreateFunction(TPZCompElDisc::CreateT2Disc);
  TPZGeoElQ2d::SetCreateFunction(TPZCompElDisc::CreateQ2Disc);
  gmesh->BuildConnectivity();
  TPZMaterial *mater = new TPZMatPoisson3d(100,3);
  mater->SetForcingFunction(Function);
  cmesh->InsertMaterialObject(mater);
  //condi��es de contorno
  //if(problemtype == 1) CC_Projecao(mater);
  //if(problemtype == 2) CC_Poisson(mater);
  int ordem=2;
  //cout << "\nEntre ordem -> ";
  //cin >> ordem;
//  TPZCompEl::gOrder = ordem;
  cmesh->SetDefaultOrder(ordem);
  cout << endl;
  cmesh->AutoBuild();
}

void Prisma(){

  CriacaoDeNos(6,prisma);
  //elemento de volume
  TPZVec<int> nodes;
  nodes.Resize(6);
  nodes[0] = 0;
  nodes[1] = 1;
  nodes[2] = 2;
  nodes[3] = 3;
  nodes[4] = 4;
  nodes[5] = 5;
  new TPZGeoElPr3d(nodes,100,*gmesh);//TPZGeoElPr3d *elg3d =
  int interfdim = 2;
  TPZCompElDisc::gInterfaceDimension = interfdim;
  TPZGeoElPr3d::SetCreateFunction(TPZCompElDisc::CreatePr3Disc);
  TPZGeoElT2d::SetCreateFunction(TPZCompElDisc::CreateT2Disc);
  TPZGeoElQ2d::SetCreateFunction(TPZCompElDisc::CreateQ2Disc);
  gmesh->BuildConnectivity();
  TPZMaterial *mater = new TPZMatPoisson3d(100,3);
  mater->SetForcingFunction(Function);
  cmesh->InsertMaterialObject(mater);
  //condi��es de contorno
  //if(problemtype == 1) CC_Projecao(mater);
  //if(problemtype == 2) CC_Poisson(mater);
  int ordem=2;
  //cout << "\nEntre ordem -> ";
  //cin >> ordem;
//  TPZCompEl::gOrder = ordem;
  cmesh->SetDefaultOrder(ordem);
  cout << endl;
  cmesh->AutoBuild();
}

void Tetraedro(){

  CriacaoDeNos(4,tetraedro);
  //elemento de volume
  TPZVec<int> nodes;
  nodes.Resize(4);
  nodes[0] = 0;
  nodes[1] = 1;
  nodes[2] = 2;
  nodes[3] = 3;
  new TPZGeoElT3d(nodes,100,*gmesh);//TPZGeoElT3d *elg3d =
  int interfdim = 2;
  TPZCompElDisc::gInterfaceDimension = interfdim;
  //construtor descont�nuo
  TPZGeoElT3d::SetCreateFunction(TPZCompElDisc::CreateT3Disc);
  TPZGeoElPi3d::SetCreateFunction(TPZCompElDisc::CreatePi3Disc);
  TPZGeoElQ2d::SetCreateFunction(TPZCompElDisc::CreateQ2Disc);
  TPZGeoElT2d::SetCreateFunction(TPZCompElDisc::CreateT2Disc);
  gmesh->BuildConnectivity();
  TPZMaterial *mater = new TPZMatPoisson3d(100,3);
  mater->SetForcingFunction(Function);
  cmesh->InsertMaterialObject(mater);
  int ordem=2;
  //cout << "\nEntre ordem -> ";
  //cin >> ordem;
//  TPZCompEl::gOrder = ordem;
  cmesh->SetDefaultOrder(ordem);
  cout << endl;
  cmesh->AutoBuild();
}

void CC_Projecao(TPZMaterial *mater){

  cout << "main nao usar\n";

   TPZBndCond *bc;
   //   REAL big = 1.e12;
   TPZFMatrix val1(1,1,0.),val2(1,1,0.);
   TPZGeoEl *elg0 = gmesh->ElementVec()[0];

   //CC na face 20
   val1.Zero();//para uso CC tipo 2
   val2.Zero();
   TPZGeoElBC(elg0,20,-1,*gmesh);
   bc = mater->CreateBC(-1,0,val1,val2);
   bc->SetForcingFunction(Function);
   cmesh->InsertMaterialObject(bc);

   //CC na face 21
   TPZGeoElBC(elg0,21,-2,*gmesh);
   bc = mater->CreateBC(-2,0,val1,val2);
   bc->SetForcingFunction(Function);
   cmesh->InsertMaterialObject(bc);

   //CC na face 22
   TPZGeoElBC(elg0,22,-3,*gmesh);
   bc = mater->CreateBC(-3,0,val1,val2);
   bc->SetForcingFunction(Function);
   cmesh->InsertMaterialObject(bc);

   //CC na face 23
   TPZGeoElBC(elg0,23,-4,*gmesh);
   bc = mater->CreateBC(-4,0,val1,val2);
   bc->SetForcingFunction(Function);
   cmesh->InsertMaterialObject(bc);

   //CC na face 24
   TPZGeoElBC(elg0,24,-5,*gmesh);
   bc = mater->CreateBC(-5,0,val1,val2);
   bc->SetForcingFunction(Function);
   cmesh->InsertMaterialObject(bc);

   //CC na face 25
   TPZGeoElBC(elg0,25,-6,*gmesh);
   bc = mater->CreateBC(-6,0,val1,val2);
   bc->SetForcingFunction(Function);
   cmesh->InsertMaterialObject(bc);
}

void CC_Poisson(TPZMaterial *mater){

  cout << "main nao usar\n";

  //nulo nos lados contidos nos eixos
  TPZBndCond *bc;
  //  REAL big = 1.e12;
  TPZFMatrix val1(1,1,0.),val2(1,1,0.);
  TPZGeoEl *elg0 = gmesh->ElementVec()[0];

  //CC na face 20
  val1.Zero();//para uso CC tipo 2
  val2.Zero();
  TPZGeoElBC(elg0,20,-1,*gmesh);
  bc = mater->CreateBC(-1,0,val1,val2);
  bc->SetForcingFunction(BCPoisson);
  cmesh->InsertMaterialObject(bc);

  //CC na face 21
  TPZGeoElBC(elg0,21,-2,*gmesh);
  bc = mater->CreateBC(-2,0,val1,val2);
  bc->SetForcingFunction(BCPoisson);
  cmesh->InsertMaterialObject(bc);

  //CC na face 22
  TPZGeoElBC(elg0,22,-3,*gmesh);
  bc = mater->CreateBC(-3,1,val1,val2);
  bc->SetForcingFunction(BCPoisson);
  cmesh->InsertMaterialObject(bc);

  //CC na face 23
  TPZGeoElBC(elg0,23,-4,*gmesh);
  bc = mater->CreateBC(-4,1,val1,val2);
  bc->SetForcingFunction(BCPoisson);
  cmesh->InsertMaterialObject(bc);

  //CC na face 24
  TPZGeoElBC(elg0,24,-5,*gmesh);
  bc = mater->CreateBC(-5,0,val1,val2);
  bc->SetForcingFunction(BCPoisson);
  cmesh->InsertMaterialObject(bc);

  //CC na face 25
  TPZGeoElBC(elg0,25,-6,*gmesh);
  bc = mater->CreateBC(-6,1,val1,val2);
  bc->SetForcingFunction(BCPoisson);
  cmesh->InsertMaterialObject(bc);
}

void BCPoisson(TPZVec<REAL> &x,TPZVec<REAL> &result){

  //solu��o u = x2+y2+z2, dom�nio cubo unit�rio
  if( fabs(x[0]-1.0) < 1.e-8 ){//face 22
    result[0] = 2.0*x[0];
    return;
  }
  if( fabs(x[1]-1.0) < 1.e-8 ){//face 23
    result[0] = 2.0*x[1];
    return;
  }
  if( fabs(x[2]-1.0) < 1.e-8 ){//face 25
    result[0] = 2.0*x[2];
    return;
  }
  //faces nos planos X=0 (face 24), Y=0 (face 21), Z=0 (facec 20)
  //no plano Xi=0 a derivada normal depende de Xi, logo � nula
  //a condi��o imposta � de Dirichlet, a CC de Neumman � nula
  if( fabs(x[0]) < 1.e-8 || fabs(x[1]) < 1.e-8 || fabs(x[2]) < 1.e-8)
    result[0] = x[0]*x[0]+x[1]*x[1]+x[2]*x[2];
}

// void Function(TPZVec<REAL> &x,TPZVec<REAL> &result){
//   //sol u = x*x + y*y + z*z , -laplaciano(u) = f
//   if(problemtype == 1) result[0]  = x[0]*x[0]+x[1]*x[1]+x[2]*x[2];
//   if(problemtype == 2) result[0] = -6.0;//Laplaciano
// }

void Solution(TPZVec<REAL> &x,TPZVec<REAL> &result,TPZFMatrix &deriv){

  Function(x,result);

  deriv(0,0) = 0.;
  deriv(1,0) = 0.;
  deriv(2,0) = 0.;

  if(problem==0) return;

  if(problem==1) {
    deriv(0,0) =  1.0;
    return;
  }

  if(problem==10){
    int i,nel = x.NElements();
    REAL norma = 0.;
    for(i=0;i<nel;i++) norma += (x[i]-x0[i])*(x[i]-x0[i]);
    norma = sqrt(norma);//norma(x-x0)
    if(norma <= r){
      REAL teta = norma*pi/r;
      REAL val = -(2.0*pi/r)*sin(teta) / norma;
      deriv(0,0) = val*(x[0]-x0[0]);
      deriv(1,0) = val*(x[1]-x0[1]);
      deriv(2,0) = val*(x[2]-x0[2]);
      return;
    }
    return;//derivada nula
  }
}

void Divisao (TPZCompMesh *cmesh,int key){

  if(key < 0) return;
  TPZVec<int> csub(0);
  int n1=1;
  while(n1) {
    cout << "\nId do elemento geometrico a dividir ? : ";
    cin >> n1;
    if(n1 < 0) break;
    int nelc = cmesh->ElementVec().NElements();
    int el=0;
    TPZCompEl *cpel=0;
    for(el=0;el<nelc;el++) {
      cpel = cmesh->ElementVec()[el];
      if(cpel && cpel->Reference()->Id() == n1) break;
    }
    if(cpel && el < nelc && cpel->Type() == 16){
      PZError << "main::Divisao elemento interface (nao foi dividido!)\n\n";
      cout << "Elementos divissiveis:\n";
      for(el=0;el<nelc;el++) {
	cpel = cmesh->ElementVec()[el];
	if(cpel && cpel->Type() != 16){
	  TPZGeoEl *gel = cpel->Reference();
	  if(gel) cout << gel->Id() << ",";
	}
      }
    } else {
      if(!el || el < nelc) cmesh->Divide(el,csub,0);
      else {
	cout << "main::Divisao elemento sem referencia\n";
	ContagemDeElementos();
      }
      n1 = 1;
    }
  }
}

void TestShapeDescontinous(){

  int nel = cmesh->ElementVec().NElements(),i;
  REAL C[3];
  if(0){
    for(i=0;i<nel;i++){
      //      TPZCompEl *comp = (cmesh->ElementVec()[i]);
      //      TPZCompElDisc *disc = dynamic_cast<TPZCompElDisc *> (comp);
      //      C[i] = disc->Constant();
    }
  }
  C[0] = 1.;
  C[1] = 1.;
  C[2] = 1.;
  TPZVec<REAL> X0(3,0.5),X(3,0.2);
  X[1] = 0.3;
  X[2] = 0.4;
  int degree = 3;
  TPZFMatrix phi,dphi;
  //  int const N=1;
  TPZShapeDisc::Shape1D(C[0],X0,X,degree,phi,dphi);
  phi.Print("Uni-dimensional",cout);
  dphi.Print("Uni-dimensional",cout);
  phi.Resize(0,0);
  dphi.Resize(0,0);
  TPZShapeDisc::Shape2D(C[1],X0,X,degree,phi,dphi);
  phi.Print("Bi-dimensional",cout);
  dphi.Print("Bi-dimensional",cout);
  phi.Resize(0,0);
  dphi.Resize(0,0);
  TPZShapeDisc::Shape3D(C[2],X0,X,degree,phi,dphi);
  phi.Print("Tri-dimensional",cout);
  dphi.Print("Tri-dimensional",cout);
}

void CycleRefinements(TPZCompMesh& cm, int numcycles, int minel, int maxel, ofstream &/*out*/){

  //  int maxBC=0;
  int numel = cm.NElements();
  if(numel > minel) minel = numel;
  TPZAdmChunkVector<TPZCompEl *> &elemvec = cm.ElementVec();
  int ic,aleat;
  cout <<  "\nEntre numero inteiro aleatorio ";// << 1;
  //aleat = 15;
  cin >> aleat;
  srand(aleat);
  clock_t start,end,begin,ttot=0;
  ofstream divide("adapta.out");
  int divididos,agrupados,totdiv=0,totagr=0,realdiv=0;
  //  int elbc=0;
  begin = clock();
  for(ic=0; ic<numcycles; ic++){
    start = clock();
    cout << "\nProcessando ciclo numero (inicio) : " << (ic+1) << endl;
    numel = cm.NElements() - cm.ElementVec().NFreeElements();
    divididos = 0;
    while(numel < maxel) {
      int elindex = rand()%cm.NElements();
      TPZCompEl *cel = elemvec[elindex];
      if(cel){
	if(cel->Type() == 16) continue;
/* 	static int key = 1; */
/* 	if(key){ */
/* 	  cout << "\nNumero maximo de elementos de contorno?  "; */
/* 	  cin >> maxBC; */
/* 	  key = 0; */
/* 	} */
/* 	if(cel->Reference()->MaterialId() < 0 && ++elbc > maxBC) continue; */
	divide << "divide id   =  " << cel->Reference()->Id();
	divide << "\telindex    =  " << elindex << endl;
	divide.flush();
	if(cel->Reference()->HasSubElement()){
	  PZError << "Elemento ja' era dividido\n";
	} else {
	  if(cel->Reference()->MaterialId() < 0){
	    PZError << "Elemento BC nao era dividido\n";
	  }
	  else PZError << "Elemento nao era dividido\n";
	  realdiv++;//todos os divididos que n�o eram agrupados
	}
	TPZVec<int> subindex;
	cel->Divide(elindex,subindex);
	divididos++;
	numel = cm.NElements() - cm.ElementVec().NFreeElements();
      }
    }
    totdiv += divididos;
    //}//TESTE
    //return;//TESTE
    //{//TESTE
    end = clock();
    cout << "\nEnd divide cycle: "  << (ic+1) << endl;
    cout << "Total divididos do ciclo " << (ic+1) << " -> " << divididos << endl;
    clock_t segundos = ((end - start)/CLOCKS_PER_SEC);
    cout << segundos << " segundos" << endl;
    cout << segundos/60.0 << " minutos" << endl << endl;
    ttot += segundos;
    start = end;
    //cm.Print(out);
    //int i;
    agrupados = 0;
    while(numel > minel) {
      int elindex = rand()%cm.NElements();
      TPZCompEl *cel = cm.ElementVec()[elindex];
      if(cel) {
	TPZVec<int> subindex;
	int ID = cel->Reference()->Id();//depois cel pode ser apagado
	if(IsGroup(elindex,cel,subindex)){
	  divide << "agrupa id   =  " << ID;
	  divide << "\tpai id   =  " << cel->Reference()->Father()->Id();
	  divide << "        elindex    =  " << elindex << endl;
	  divide.flush();
	  //out << "Coarsening "; for(i=0;i<4;i++) out << subindex[i] << "/"
	  //    << cm.ElementVec()[subindex[i]]->Reference()->Id() << " ";
	  //out << endl;
	  cm.Coarsen(subindex,elindex);
	  agrupados++;
	  //out << "Created " << elindex << "/"
	  //<< cm.ElementVec()[elindex]->Reference()->Id() << endl;
	  cout << "Agrupado " << ID << endl;
	  //out.flush();
	  cm.CleanUpUnconnectedNodes();
	  numel = cm.NElements()-cm.ElementVec().NFreeElements();
	}
      }
    }//
    totagr += agrupados;
    cout << "\nEnd coarsen cycle: " << (ic+1) << endl;
    cout << "Total agrupados do ciclo " << (ic+1) << " -> " <<  agrupados << endl;
    cout << "\nCleanUpUnconnectedNodes()\n";
    cm.CleanUpUnconnectedNodes();
    end = clock();
    segundos = ((end - start)/CLOCKS_PER_SEC);
    cout << segundos << " segundos" << endl;
    cout << segundos/60.0 << " minutos" << endl << endl;
    ttot += segundos;
    //out.flush();
    //cm.Print(out);
    cout << "\nTotal divididos = " << totdiv;
    cout << "\nTotal realmente divididos = " << realdiv;
    cout << "\nTotal agrupados = " << totagr << endl;
    end = clock();
    cout << "\nEnd cycle : total time: " << endl;
    segundos = ((end - begin)/CLOCKS_PER_SEC);
    cout << segundos << " segundos" << endl;
    cout << segundos/60.0 << " minutos" << endl;
    cout << ttot << " segundos (soma das parciais)" << endl;
    cout << ttot/60.0 << " minutos (soma das parciais)" << endl << endl;
  }
  divide << "\n\n\t\tFIM CICLO DE REFINAMENTOS E AGRUPAMENTOS\n\n\n";
  divide.close();
}

int IsGroup(int /*ind*/, TPZCompEl *cel, TPZVec<int> &indgroup) {

  TPZGeoEl *gel = cel->Reference();
  if(!gel) return 0;
  TPZGeoEl *gelf = gel->Father();
  if(!gelf) return 0;
  int nsub = gelf->NSubElements();
  TPZVec<TPZGeoEl *> subel(nsub);
  int is;
  for(is=0; is<nsub; is++) subel[is] = gelf->SubElement(is);
  TPZVec<TPZCompEl *> csubel(nsub);
  indgroup.Resize(nsub);
  for(is=0; is<nsub; is++) {
    csubel[is] = subel[is]->Reference();
    if(csubel[is] == 0) return 0;
    indgroup[is] = csubel[is]->Index();
  }
  return 1;
}

void ContagemDeElementos(){

  int poin=0,line=0,tria=0,quad=0,tetr=0,pira=0,pris=0,hexa=0,disc=0,inte=0;
  int nelem = cmesh->ElementVec().NElements();
  int k,totel=0,bcel=0,niv = 0,nivmax=0;
  for(k=0;k<nelem;k++){
    TPZCompEl *comp = cmesh->ElementVec()[k];
    if(!comp) continue;
    totel++;
    if(comp->Reference()->MaterialId() < 0) bcel++;
    niv = comp->Reference()->Level();
    if(nivmax < niv) nivmax = niv;
/*     if(comp->Type() == 00) poin++; */
/*     if(comp->Type() == 01) line++; */
/*     if(comp->Type() == 02) tria++; */
/*     if(comp->Type() == 03) quad++; */
/*     if(comp->Type() == 04) tetr++; */
/*     if(comp->Type() == 05) pira++; */
/*     if(comp->Type() == 06) pris++; */
/*     if(comp->Type() == 07) hexa++; */
    if(comp->Type() == 15) disc++;
    if(comp->Type() == 16) inte++;
  }
  nelem = gmesh->ElementVec().NElements();
  int total=0,nivmax2=0;
  for(k=0;k<nelem;k++){
    TPZGeoEl *geo = gmesh->ElementVec()[k];
    if(!geo) continue;
    total++;
    niv = geo->Level();
    if(nivmax2 < niv) nivmax2 = niv;
    if(geo->Reference()){
      int nsides = geo->NSides();
      if(nsides ==  1) poin++;
      if(nsides ==  3) line++;
      if(nsides ==  7) tria++;
      if(nsides ==  9) quad++;
      if(nsides == 15) tetr++;
      if(nsides == 19) pira++;
      if(nsides == 21) pris++;
      if(nsides == 27) hexa++;
    }
  }
  REAL passo = dynamic_cast<TPZConservationLaw *>(mat)->TimeStep();
  cout << "\nTotal de elementos computacionais  : " << totel;
  cout << "\nTotal de elementos de dominio      : " << abs(totel-bcel);
  cout << "\nTotal de elementos de contorno     : " << bcel;
  cout << "\nTotal de elementos ponto           : " << poin;
  cout << "\nTotal de elementos linha           : " << line;
  cout << "\nTotal de elementos triangulo       : " << tria;
  cout << "\nTotal de elementos quadrilatero    : " << quad;
  cout << "\nTotal de elementos tetraedro       : " << tetr;
  cout << "\nTotal de elementos piramide        : " << pira;
  cout << "\nTotal de elementos prisma          : " << pris;
  cout << "\nTotal de elementos hexaedro        : " << hexa;
  cout << "\nTotal do tipo discontinuo          : " << disc;
  cout << "\nTotal de tipo interface            : " << inte;
  cout << "\nTotal de nos                       : " << gmesh->NodeVec().NElements();
  cout << "\nTotal de elementos geometricos     : " << total;
  cout << "\nTamanho do vetor de connects       : " << cmesh->NConnects();
  cout << "\nTamanho do vetor el. comput.       : " << cmesh->NElements();
  cout << "\nNivel maximo comput. atingido      : " << nivmax;
  cout << "\nGrau de interpola��o               : " << grau;
  cout << "\nCFL                                : " << CFL;
  cout << "\nPasso de tempo                     : " << passo;
  cout << "\nNivel maximo geomet. atingido      : " << nivmax << endl << endl;
}

int Nivel(TPZGeoEl *gel);
void NivelDivide(TPZCompMesh *cmesh){

  TPZVec<int> csub(0);
  //int nivel;
  cout << "\nmain::Divisao todos os elementos da malha serao divididos!\n";
  //cout << "\nmain::Divisao Nivel da malha final ? : ";
  //cin >> nivel;
  cout << "\nNivel da malha a ser atingido = " << nivel << endl;
  int nelc = cmesh->ElementVec().NElements();
  int el,actual;
  TPZCompEl *cpel;
  TPZGeoEl *gel;
  el = -1;
  while(++el<nelc) {
    cpel = cmesh->ElementVec()[el];
    if(!cpel) continue;
    if(cpel->Type() == 16) continue;
    if(cpel->Material()->Id() < 0) continue;
    gel = cpel->Reference();
    actual = Nivel(gel);
    if(actual < nivel){
      cmesh->Divide(el,csub,0);
      nelc = cmesh->ElementVec().NElements();
      el = -1;
      continue;
    }
  }
}

int Nivel(TPZGeoEl *gel){
  //retorna o n�vel do elemento gel
  if(!gel) return -1;
  TPZGeoEl *fat = gel->Father();
  if(!fat) return 0;
  int niv = 0;
  while(fat){
    fat = fat->Father();
    niv++;
  }
  return niv;
}

static REAL point[5][3] = {{-.8,0.,0.},{-.4,.0,.0},{.0,.0,.0},{.4,.0,.0},{.8,0.,0.}};//linha
static REAL quad[5][3] = {{-.5,-.5,0.},{.5,-.5,.0},{.5,.5,.0},{-.5,.5,.0},{.0,.0,0.}};//quadrilatero
static REAL hexa[9][3] = { {-0.8,-0.8,-0.8},{0.8,-0.8,-0.8},{0.8,0.8,-0.8},{-0.8,0.8,-0.8},
			   {-0.8,-0.8,00.8},{0.8,-0.8,00.8},{0.8,0.8,00.8},{-0.8,0.8,00.8},{0.,0.,0.} };//hexaedro
void PostProcess(TPZGeoMesh &gmesh,ostream &out) {

  int nel = gmesh.Reference()->ElementVec().NElements();
  if(nel > 1000){
    cout << "main::PostProcess mas de 10000 elementos -> processa 2000\n";
  }
  int idmax = 0,dim,chega=1,finish=-1;
  for(int iel=0;iel<nel;iel++){//procurando o id mais alto da lista
    if(++finish >= 2000) return;
    TPZCompEl *cel = gmesh.Reference()->ElementVec()[iel];
    if(!cel) continue;
    TPZGeoEl *el = gmesh.ElementVec()[iel];
    if(chega && cel->Type()==15) {dim = el->Dimension(); chega = 0;}
    int id = el->Id();
    if(id > idmax) idmax = id;
  }
  for(int iel=0;iel<nel;iel++) {
    if(!gmesh.Reference()->ElementVec()[iel]) continue;
    int elemtype = gmesh.Reference()->ElementVec()[iel]->Type();
    if(elemtype==16) continue;
    TPZCompEl *el = gmesh.Reference()->ElementVec()[iel];
    if(el->Material()->Id() < 0) continue;
    TPZGeoEl *gel = el->Reference();
    if(el && gel) {
      TPZGeoElPoint  *el0d=0;
      TPZGeoEl1d  *el1d=0;
      TPZGeoElT2d *elt2d=0;
      TPZGeoElQ2d *elq2d=0;
      TPZGeoElT3d *elt3d=0;
      TPZGeoElPi3d *elpi3d=0;
      TPZGeoElPr3d *elpr3d=0;
      TPZGeoElC3d  *elc3d=0;
      if(elemtype==0) el0d   = (TPZGeoElPoint *) gel;
      if(elemtype==1) el1d   = (TPZGeoEl1d    *) gel;
      if(elemtype==2) elt2d  = (TPZGeoElT2d   *) gel;
      if(elemtype==3) elq2d  = (TPZGeoElQ2d   *) gel;
      if(elemtype==4) elt3d  = (TPZGeoElT3d   *) gel;
      if(elemtype==5) elpi3d = (TPZGeoElPi3d  *) gel;
      if(elemtype==6) elpr3d = (TPZGeoElPr3d  *) gel;
      if(elemtype==7) elc3d  = (TPZGeoElC3d   *) gel;
      int nsides = gel->NSides();
      if(elemtype==15){
	if(nsides==1) el0d   = (TPZGeoElPoint *) gel;
	if(nsides==3) el1d   = (TPZGeoEl1d    *) gel;
	if(nsides==7) elt2d  = (TPZGeoElT2d   *) gel;
	if(nsides==9) elq2d  = (TPZGeoElQ2d   *) gel;
	if(nsides==27) elc3d  = (TPZGeoElC3d   *) gel;
      }
      out << "Elemento " << el->Reference()->Id() << endl;;
      TPZManVector<REAL> sol(1);
      TPZVec<REAL> csi(3,0.),x(3);
      int np = 5;
      if(dim==3) np = 9;
      for(int p=0;p<np;p++) {
	if(dim==1){
	  csi[0] = point[p][0];
	  csi[1] = point[p][1];
	  csi[2] = point[p][2];
	}
	if(dim==2){
	  csi[0] = quad[p][0];
	  csi[1] = quad[p][1];
	  csi[2] = quad[p][2];
	}
	if(dim==3){
	  csi[0] = hexa[p][0];
	  csi[1] = hexa[p][1];
	  csi[2] = hexa[p][2];
	}
	gel->X(csi,x);
	if(elemtype==0) el0d->Reference()->Solution(csi,0,sol);
	if(elemtype==1) el1d->Reference()->Solution(csi,0,sol);
	if(elemtype==2) elt2d->Reference()->Solution(csi,0,sol);
	if(elemtype==3) elq2d->Reference()->Solution(csi,0,sol);
	if(elemtype==4) elt3d->Reference()->Solution(csi,0,sol);
	if(elemtype==5) elpi3d->Reference()->Solution(csi,0,sol);
	if(elemtype==6) elpr3d->Reference()->Solution(csi,0,sol);
	if(elemtype==7) elc3d->Reference()->Solution(csi,0,sol);
	if(elemtype==15){
	  if(nsides==1) el0d->Reference()->Solution(csi,0,sol);
	  if(nsides==3) el1d->Reference()->Solution(csi,0,sol);
	  if(nsides==7) elt2d->Reference()->Solution(csi,0,sol);
	  if(nsides==9) elq2d->Reference()->Solution(csi,0,sol);
	  if(nsides==27) elc3d->Reference()->Solution(csi,0,sol);
	}
	out << "solucao em x    = " << x[0] << ' ' << x[1] << ' ' << x[2] << endl;
	out << "               u = " << sol[0] << endl;
      }
    }
  }
}

void Ordena(TPZVec<REAL> &coordx,TPZVec<int> &sort);
void FileNB(TPZGeoMesh &gmesh,ostream &out) {

  //  int numpoints;
  int nel = gmesh.Reference()->ElementVec().NElements();
  int idmax = 0,dim,chega=1,finish=-1;
  for(int iel=0;iel<nel;iel++){//procurando o id mais alto da lista
    if(++finish >= 5000) return;
    TPZCompEl *cel = gmesh.Reference()->ElementVec()[iel];
    if(!cel) continue;
    TPZGeoEl *el = gmesh.ElementVec()[iel];
    if(chega && cel->Type()==15) {dim = el->Dimension(); chega = 0;}
    int id = el->Id();
    if(id > idmax) idmax = id;
  }
  int cap = nel*4;
  TPZVec<REAL> coordx(cap,0.),coordy(cap,0.),coordz(cap,0.);
  int count = -1,capacity=0;
  while(count++<idmax){
    for(int iel=0;iel<nel;iel++) {
      if(!gmesh.Reference()->ElementVec()[iel]) continue;
      int elemtype = gmesh.Reference()->ElementVec()[iel]->Type();
      if(elemtype==16) continue;//interface
      TPZCompEl *el = gmesh.Reference()->ElementVec()[iel];
      if(el->Material()->Id() < 0) continue;
      //s� elementos de volume
      TPZGeoEl *gel = el->Reference();
      if(el && gel) {
	if(gel->Id()==count){
	  TPZGeoElPoint  *el0d=0;
	  TPZGeoEl1d  *el1d=0;
	  TPZGeoElT2d *elt2d=0;
	  TPZGeoElQ2d *elq2d=0;
	  TPZGeoElT3d *elt3d=0;
	  TPZGeoElPi3d *elpi3d=0;
	  TPZGeoElPr3d *elpr3d=0;
	  TPZGeoElC3d  *elc3d=0;
	  if(elemtype==0) el0d   = (TPZGeoElPoint    *) gel;
	  if(elemtype==1) el1d   = (TPZGeoEl1d    *) gel;
	  if(elemtype==2) elt2d  = (TPZGeoElT2d   *) gel;
	  if(elemtype==3) elq2d  = (TPZGeoElQ2d   *) gel;
	  if(elemtype==4) elt3d  = (TPZGeoElT3d   *) gel;
	  if(elemtype==5) elpi3d = (TPZGeoElPi3d  *) gel;
	  if(elemtype==6) elpr3d = (TPZGeoElPr3d  *) gel;
	  if(elemtype==7) elc3d  = (TPZGeoElC3d   *) gel;
	  int nsides = gel->NSides();
	  if(elemtype==15){
	    if(nsides==1) el0d   = (TPZGeoElPoint *) gel;
	    if(nsides==3) el1d   = (TPZGeoEl1d    *) gel;
	    if(nsides==7) elt2d  = (TPZGeoElT2d   *) gel;
	    if(nsides==9) elq2d  = (TPZGeoElQ2d   *) gel;
	  }
	  TPZManVector<REAL> sol(1);
	  TPZVec<REAL> csi(3,0.),x(3);
	  for(int p=0;p<4;p++) {
	    if(dim==1){
	      csi[0] = point[p][0];
	      csi[1] = point[p][1];
	      csi[2] = point[p][2];
	    }
	    if(dim==2){
	      csi[0] = quad[p][0];
	      csi[1] = quad[p][1];
	      csi[2] = quad[p][2];
	    }
	    gel->X(csi,x);
	    if(elemtype==0) el0d->Reference()->Solution(csi,0,sol);
	    if(elemtype==1) el1d->Reference()->Solution(csi,0,sol);
	    if(elemtype==2) elt2d->Reference()->Solution(csi,0,sol);
	    if(elemtype==3) elq2d->Reference()->Solution(csi,0,sol);
	    if(elemtype==4) elt3d->Reference()->Solution(csi,0,sol);
	    if(elemtype==5) elpi3d->Reference()->Solution(csi,0,sol);
	    if(elemtype==6) elpr3d->Reference()->Solution(csi,0,sol);
	    if(elemtype==7) elc3d->Reference()->Solution(csi,0,sol);
	    if(elemtype==15){
	      if(nsides==1) el0d->Reference()->Solution(csi,0,sol);
	      if(nsides==3) el1d->Reference()->Solution(csi,0,sol);
	      if(nsides==7) elt2d->Reference()->Solution(csi,0,sol);
	      if(nsides==9) elq2d->Reference()->Solution(csi,0,sol);
	    }
	    if(dim==1){
	      coordx[capacity] = x[0];
	      coordy[capacity] = sol[0];
	    }
	    if(dim==2){
	      coordx[capacity] = x[0];
	      coordy[capacity] = x[1];
	      coordz[capacity] = sol[0];
	  }
	    capacity++;
	  }
	  //out.close();
	} else {
	  continue;
	}
      }
    }
  }
  if(dim==1){
    //out << "<<Graphics`MultipleListPlot`\n";
    out << "GRAPH = {";
    int k,linha = 0;
    TPZVec<int> sort(capacity);
    Ordena(coordx,sort);
    for(k=0;k<(capacity-1);k++){
      out <<  "{" << coordx[k] << "," <<  coordy[sort[k]] << "},";
      linha++;
      if(linha == 8){
	out << endl;
	linha = 0;
      }
    }
    out <<  "{" << coordx[k] << "," <<  coordy[sort[k]] << "}};";
    out << "\nListPlot[GRAPH,PlotJoined->True]";
  }
  if(dim==2){
    //out << "<<Graphics`MultipleListPlot`\n";
    out << "GRAPH = {";
    int k,linha = 0;
    TPZVec<int> sort1(capacity);//,sort2(capacity);
    Ordena(coordx,sort1);
//     TPZVec<REAL> sort3(coordy),sort4(coordz);
//     for(k=0;k<capacity;k++){
//       coordy[k] = sort3[sort1[k]];
//       coordz[k] = sort4[sort1[k]];
//     }
//    Ordena(coordy,sort2);
    for(k=0;k<(capacity-1);k++){
      out <<  "{" << coordx[k] << "," << coordy[sort1[k]] << "," << coordz[sort1[k]] << "},";
      linha++;
      if(linha == 6){
	out << endl;
	linha = 0;
      }
    }
    out <<  "{" << coordx[k] << "," << coordy[sort1[k]] << "," << coordz[sort1[k]] << "}};";
    out << "\nListPlot3D[GRAPH]";
  }
}

void Ordena(TPZVec<REAL> &coordx,TPZVec<int> &sort){

  int i,j,cap=sort.NElements();
  for(i=0;i<cap;i++) sort[i] = i;
  for(i=0;i<cap;i++){
    REAL x = coordx[i];
    for(j=i+1;j<cap;j++){
      if(coordx[j] < x){
	coordx[i] = coordx[j];
	coordx[j] = x;
	x = coordx[i];
	int aux = sort[i];
	sort[i] = sort[j];
	sort[j] = aux;
      }
    }
  }
}

void CoutTime(clock_t &start){
    end = clock();
    cout << "\nFim da etapa : "  <<  endl;
    clock_t segundos = ((end - start)/CLOCKS_PER_SEC);
    cout << segundos << " segundos" << endl;
    cout << segundos/60.0 << " minutos" << endl << endl;
}

TPZMaterial *HexaedroNovo(int ordem){

  CriacaoDeNos(8,hexaedro);
  //elemento de volume
  TPZVec<int> nodes;
  nodes.Resize(8);
  nodes[0] = 0;
  nodes[1] = 1;
  nodes[2] = 2;
  nodes[3] = 3;
  nodes[4] = 4;
  nodes[5] = 5;
  nodes[6] = 6;
  nodes[7] = 7;
  TPZGeoElC3d *elgc3d = new TPZGeoElC3d(nodes,1,*gmesh);
  //construtor descont�nuo
  TPZGeoElC3d::SetCreateFunction(TPZCompElDisc::CreateC3Disc);//de volume
  TPZGeoElQ2d::SetCreateFunction(TPZCompElDisc::CreateQ2Disc);//n�o de interface
  int interfdim = 2;
  TPZCompElDisc::gInterfaceDimension = interfdim;
  gmesh->BuildConnectivity();
  TPZVec<REAL> B(3);
  int test = 0;
  if(!problem){
    B[0] = 1.0;
    B[1] = 0.0;
    B[2] = 0.0;
  }
  if(problem==1){
    B[0] = 1.0;
    B[1] = 0.0;
    B[2] = 0.0;
  }
  if(problem==9){
    B[0] = 1.0;
    B[1] = 0.0;
    B[2] = 0.0;
    test = 0;
  }
  if(problem==10){
    test = 2;
  }
  int nummat = 1;
  int dif_artif = 1;//SUPG
  cout << "\nmain::Divisao Nivel final da malha ? : ";
  cin >> nivel;
  REAL cfl = ( 1./(2.0*(REAL)grau+1.0) );
  REAL delta_x = ( 1.0 / pow(2.0,(REAL)nivel) );//tamanho do elemento
  REAL delta = ( (10./3.)*cfl*cfl - (2./3.)*cfl + 1./10. );
  REAL maxflinha = 1.0;
  REAL delta_t = cfl*delta_x / maxflinha;//delta_t � <= que este valor
  cout << "\nMaximo de f' = " << maxflinha
       << "\nCFL = " << cfl
       << "\nDominio [0,1]x[0,1]x[0,1]"
       << "\nDelta x = " << delta_x
       << "\ndelta t = " << delta_t
       << "\ndelta = " << delta << endl;
  int dim = 3;
  TPZMaterial *mater = new TPZConsLawTest(nummat,B,dif_artif,delta_t,dim,delta,test);
  //((TPZConservationLaw *)mat)->SetIntegDegree(grau);
  mater->SetForcingFunction(Function);
  cmesh->InsertMaterialObject(mat);
  //condi��es de contorno
  TPZBndCond *bc;
  //  REAL big = 1.e12;
  TPZFMatrix val1(1,1,0.),val2(1,1,0.);

  //a �nica fun��o atual destes elemento
  //� a cria��o dos elementos de interface
  //CC na aresta 4
  val1.Zero();//para uso CC tipo 2
  val2.Zero();

  if(problem==0) val2(0,0) = 1.;
  //CC na aresta 20
  TPZGeoElBC(elgc3d,20,-1,*gmesh);
  bc = mater->CreateBC(-1,1,val1,val2);
  //bc->SetForcingFunction(Function);
  cmesh->InsertMaterialObject(bc);

  //CC na aresta 21
  TPZGeoElBC(elgc3d,21,-2,*gmesh);
  bc = mater->CreateBC(-2,1,val1,val2);
  //bc->SetForcingFunction(Function);
  cmesh->InsertMaterialObject(bc);

  //CC na aresta 22
  TPZGeoElBC(elgc3d,22,-3,*gmesh);
  bc = mater->CreateBC(-3,1,val1,val2);
  //bc->SetForcingFunction(Function);
  cmesh->InsertMaterialObject(bc);

  //CC na aresta 23
  TPZGeoElBC(elgc3d,23,-4,*gmesh);
  bc = mater->CreateBC(-4,1,val1,val2);
  //bc->SetForcingFunction(Function);
  cmesh->InsertMaterialObject(bc);

  //CC na aresta 24
  TPZGeoElBC(elgc3d,24,-5,*gmesh);
  bc = mater->CreateBC(-5,1,val1,val2);
  //bc->SetForcingFunction(Function);
  cmesh->InsertMaterialObject(bc);

  //CC na aresta 25
  TPZGeoElBC(elgc3d,25,-6,*gmesh);
  bc = mater->CreateBC(-6,1,val1,val2);
  //bc->SetForcingFunction(Function);
  cmesh->InsertMaterialObject(bc);

//   int ordem=2;
//   cout << "\nEntre ordem -> ";
//   //cin >> ordem;
//   TPZCompEl::gOrder = ordem;
  cout << endl;
  cmesh->AutoBuild();

  return mater;
}

TPZMaterial *QuadrilateroNovo(int ordem){

  CriacaoDeNos(4,quadrilatero);
  //elemento de volume
  TPZVec<int> nodes;
  nodes.Resize(4);
  nodes[0] = 0;
  nodes[1] = 1;
  nodes[2] = 2;
  nodes[3] = 3;
  TPZGeoElQ2d *elgq2d = new TPZGeoElQ2d(nodes,1,*gmesh);
  TPZGeoElQ2d::SetCreateFunction(TPZCompElDisc::CreateQ2Disc);//de volume
  TPZGeoEl1d::SetCreateFunction(TPZCompElDisc::Create1dDisc);//n�o de interface
  int interfdim = 1;
  TPZCompElDisc::gInterfaceDimension = interfdim;
  gmesh->BuildConnectivity();
  TPZVec<REAL> B(2);
  int test = 0;
  if(!problem){
    B[0] = 1.0;//1d e 2d
    B[1] = 0.0;//test=0
  }
  if(problem==1){
    B[0] = 1.0;//1d e 2d
    B[1] = 0.0;//test=0
  }
  if(problem==9){
    B[0] = 1.0;//1d e 2d
    B[1] = 0.0;//test=0
    test = 0;
  }
  if(problem==10){
    test = 1;
  }
  int nummat = 1;
  int dif_artif = 1;//SUPG
  cout << "\nmain::Divisao Nivel final da malha ? : ";
  cin >> nivel;
  REAL cfl = ( 1./(2.0*(REAL)grau+1.0) );
  REAL delta_x = ( 1.0 / pow(2.0,(REAL)nivel) );//tamanho do elemento
  REAL delta = ( (10./3.)*cfl*cfl - (2./3.)*cfl + 1./10. );
  REAL maxflinha = 1.0;
  REAL delta_t = cfl*delta_x / maxflinha;//delta_t � <= que este valor
  cout << "\nMaximo de f' = " << maxflinha
       << "\nCFL = " << cfl
       << "\nDominio [0,1]x[0,1]"
       << "\nDelta x = " << delta_x
       << "\ndelta t = " << delta_t
       << "\ndelta = " << delta << endl;
  int dim = 2;
  TPZMaterial *mater = new TPZConsLawTest(nummat,B,dif_artif,delta_t,dim,delta,test);
  //((TPZConservationLaw *)mat)->SetIntegDegree(grau);
  mater->SetForcingFunction(Function);
  cmesh->InsertMaterialObject(mat);
  //condi��es de contorno
  TPZBndCond *bc;
  //  REAL big = 1.e12;
  TPZFMatrix val1(1,1,0.),val2(1,1,0.);

  //a �nica fun��o atual destes elemento
  //� a cria��o dos elementos de interface
  //CC na aresta 4
  val1.Zero();//para uso CC tipo 2
  val2.Zero();

  if(problem==0) val2(0,0) = 1.;
  TPZGeoElBC(elgq2d,4,-1,*gmesh);
  bc = mater->CreateBC(-1,1,val1,val2);
  //bc->SetForcingFunction(Function);
  cmesh->InsertMaterialObject(bc);

  //CC na aresta 5
  TPZGeoElBC(elgq2d,5,-2,*gmesh);
  bc = mater->CreateBC(-2,1,val1,val2);
  //bc->SetForcingFunction(Function);
  cmesh->InsertMaterialObject(bc);

  //CC na aresta 6
  TPZGeoElBC(elgq2d,6,-3,*gmesh);
  bc = mater->CreateBC(-3,1,val1,val2);
  //bc->SetForcingFunction(Function);
  cmesh->InsertMaterialObject(bc);

  //CC na aresta 7
  TPZGeoElBC(elgq2d,7,-4,*gmesh);
  bc = mater->CreateBC(-4,1,val1,val2);
  //bc->SetForcingFunction(Function);
  cmesh->InsertMaterialObject(bc);

  cout << endl;
  cmesh->AutoBuild();

  return mater;
}
//,ur(1,6),dux2(6,1),duy2(6,1)
// 	 ur(0,k) = u(k,0);
// 	 dux(0,k) = du(0,k);
// 	 duy(0,k) = du(1,k);
//        Fun =  (dphi(0,i)*(dux*Kn1n1)+
// 	       dphi(0,i)*(duy*Kn1n2)+
// 	       dphi(1,i)*(dux*Kn2n1)+
// 	       dphi(1,i)*(duy*Kn2n2)+
// 	       dphi(0,i)*(ur*  Bn10 )+
// 	       dphi(1,i)*(ur*  Bn20 )+
// 	       phi(i,0) *(dux*B0n1 )+
// 	       phi(i,0) *(duy*B0n2 )+
// 	       phi(i,0) *(ur*  B000 ));
//       for(idf=0; idf<6; idf++) ef(6*i+idf,0) += weight*Fun(0,idf)*phi(i,0);
