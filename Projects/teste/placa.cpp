//$Id: placa.cpp,v 1.7 2009-06-24 20:14:55 phil Exp $
#include "TPZGeoCube.h"
#include "pzshapecube.h"
#include "TPZRefCube.h"
#include "pzshapelinear.h"
#include "TPZGeoLinear.h"
#include "TPZRefLinear.h"
#include "pzrefquad.h"
#include "pzshapequad.h"
#include "pzgeoquad.h"
#include "pzshapetriang.h"
#include "pzreftriangle.h"
#include "pzgeotriangle.h"
#include "pzshapeprism.h"
#include "pzrefprism.h"
#include "pzgeoprism.h"
#include "pzshapetetra.h"
#include "pzreftetrahedra.h"
#include "pzgeotetrahedra.h"
#include "pzshapepiram.h"
#include "pzrefpyram.h"
#include "pzgeopyramid.h"
#include "pzrefpoint.h"
#include "pzgeopoint.h"
#include "pzplaca.h"


#include "TPZGeoElement.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzfmatrix.h"
#include "pzbndcond.h"
//#include "pztempmat.h"
#include "pzcompel.h"
#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "pzskylmat.h"
#include "pzstepsolver.h"
#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzstack.h"
#include "pzvec.h"
#include "pzsolve.h"
//#include "pzelgpoint.h"
#include "pzelmat.h"
#include "pzelasmat.h"
#include "pzmattest.h"
#include "pzmat1dlin.h"
#include "pzmat2dlin.h"
#include "pzpoisson3d.h"
#include "pzmaterial.h"
#include "TPZCompElDisc.h"
#include "TPZShapeDisc.h"
#include "TPZInterfaceEl.h"
#include "pzreal.h"
#include "pzdxmesh.h"

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <ostream>
#include <string.h>
#include <math.h>
using namespace std;

//#include "BCTest.cc"

//static double quadrilatero[4][3]  = { {0.,0.,0.},{1.,0.,0.},{0.,1.,0.},{1.,1.,0.} };

static double NosQuads[9][3]  = { {-1.,-1.,0.},{0.,-1.,0.},{1.,-1.,0.},{-1.,0.,0.},{0.,0.,0.},
				  {1.,0.,0.},{-1.,1.,0.},{0.,1.,0.},{1.,1.,0.} };
static int    IncidQuads[4][8] = {{0,1,4,3},{1,2,5,4},{4,5,8,7},{3,4,7,6}};

static double Teste1Nos[4][3]  = { {0.,0.,0.},{2.,0.,0.},{2.,2.,0.},{0.,2.,0.}};
static int    Teste1Elems[1][8] = {{0,1,2,3}};

static double PolNos[4][3]  = { {0.,0.,0.},{1.,0.,0.},{1.,1.,0.},{0.,1.,0.}};
static int    PolElems[1][8] = {{0,1,2,3}};

static double Teste2Nos[5][3]  = { {0.,0.,0.},{2.,0.,0.},{2.,2.,0.},{0.,2.,0.},{1.,1.,0.}};
static int    Teste2Elems[4][8] = {{0,1,4},{1,2,4},{2,3,4},{3,0,4}};


static double Nodos3d[18][3] = {{-1.,-1.,0.},{0.,-1.,0.},{1.,-1.,0.},{-1.,0.,0.,},{0.,0.,0.},
				{1.,0.,0.},{-1.,1.,0.},{0.,1.,0.},{1.,1.,0.},
				{-1.,-1.,.2},{0.,-1.,.2},{1.,-1.,.2},{-1.,0.,.2},{0.,0.,.2},
				{1.,0.,.2},{-1.,1.,.2},{0.,1.,.2},{1.,1.,.2}};

static int  Cubic3d[4][8]  = { {0,1,4,3,9,10,13,12},{1,2,5,4,10,11,14,13},
			       {4,5,8,7,13,14,17,16},{3,4,7,6,12,13,16,15} };

void ProcessamentoLocal(TPZGeoMesh &gmesh,std::ostream &out);
void CriaNos(int num, TPZGeoMesh &geomesh, double list [20][3]);
void CriaElementos(int numelem,int ncon,TPZGeoMesh &geomesh, int list[20][8] );
void CriaCondCont(TPZGeoMesh &gmesh);
void CriaCondContTeste1(TPZGeoMesh &gmesh);
void CriaCondContTeste2(TPZGeoMesh &gmesh);
void CriaCondContTeste3(TPZGeoMesh &gmesh);
void CriaCondContTeste4(TPZGeoMesh &gmesh);
void ExecutaAnalysis(TPZCompMesh &cmesh,TPZMaterial *mat);
void PosProcessamento(TPZAnalysis &an);
int  Nivel(TPZGeoEl *gel);
void NivelDivide(TPZCompMesh *cmesh);
void LoadSolution(TPZFMatrix &axes,TPZVec<REAL> &X,TPZFMatrix &u,TPZFMatrix &du);
void LoadSolution1(TPZFMatrix &axes,TPZVec<REAL> &x,TPZFMatrix &u,TPZFMatrix &du);
void LoadSolution2(TPZFMatrix &axes,TPZVec<REAL> &X,TPZFMatrix &u,TPZFMatrix &du);
void LoadSolution3(TPZFMatrix &axes,TPZVec<REAL> &X,TPZFMatrix &u,TPZFMatrix &du);
void LoadSolution10(TPZFMatrix &axes,TPZVec<REAL> &x,TPZFMatrix &u,TPZFMatrix &du);
void LoadSolution11(TPZFMatrix &axes,TPZVec<REAL> &x,TPZFMatrix &u,TPZFMatrix &du);
void LoadSolution12(TPZFMatrix &axes,TPZVec<REAL> &X,TPZFMatrix &u,TPZFMatrix &du);
void LoadSolution13(TPZFMatrix &axes,TPZVec<REAL> &X,TPZFMatrix &u,TPZFMatrix &du);
void LoadSolution14(TPZFMatrix &axes,TPZVec<REAL> &X,TPZFMatrix &u,TPZFMatrix &du);
void BCSolution(TPZVec<REAL> &x,TPZVec<REAL> &result);
void Solution(TPZVec<REAL> &x,TPZVec<REAL> &result,TPZFMatrix &deriv);
REAL SolutionError(TPZCompMesh *cmesh,char *title);
REAL NormSolExact(TPZCompMesh *cmesh,char *title);
TPZMaterial *LerMaterial(char *filename,TPZCompMesh &cmesh);
void AdaptativeProcedure(REAL erro_adm,int numiter,int resolution,TPZMaterial *mat,TPZAnalysis &an);
int AdaptiveMesh(TPZCompMesh &cmesh,REAL tol,int iter);
int NumElComp(TPZCompMesh &cmesh);
ofstream erros("ERROS.out");
static int problema;
static REAL NormSol=-1.0;

/**
 * -> PROGRAMA PRINCIPAL <- -> PROGRAMA PRINCIPAL <- -> PROGRAMA PRINCIPAL <- -> PROGRAMA PRINCIPAL <-
 */

int main() {

//   TPZFMatrix axes(3,3,0.);
//   TPZVec<REAL> X(3);
//   TPZFMatrix u(6,1),du(2,6);
//   axes(0,0) = 1.0;
//   axes(1,1) = 1.0;
//   axes(2,2) = 1.0;
//   X[0] = .123;
//   X[1] = .456;
//   LoadSolution(axes,X,u,du);
//   u.Print("SOLUCAO U");
//   du.Print("SOLUCAO DU");
//   return 0;

  cout << "main::Problema : "
       << "[1:xy]\n"
       << "[12:x*x]\n"
       << "[13:y*y]\n"
       << "[14:x*x*y]\n"
       << "[2:x^p*y^q*(x-1)^p*(y-1)^q]\n"
       << "[3:exp(-(x-1)^2)*exp(-(y-1)^2)]\n"
       << "[11:x*x*y*y]\n";

  cin >> problema;
  //inserindo a ordem de interpola��o dos elementos e do espa�o
  int ord;
  cout << "Entre ordem : 1,2,3,4,5 : -> ";
  cin >> ord;
  //ord = 4; cout << endl;
//  TPZCompEl::gOrder = ord;
  TPZCompEl::SetgOrder(ord);

  //malha geometrica
  TPZGeoMesh *geomesh = new TPZGeoMesh;

  geomesh->Reference();

  //cria n�s geom�tricos
  if(0) CriaNos(9, *geomesh, NosQuads);
  if(0) CriaNos(18, *geomesh, Nodos3d);
  if(0) CriaNos(4, *geomesh, Teste1Nos);
  if(0) CriaNos(5, *geomesh, Teste2Nos);
  if(1) CriaNos(4, *geomesh, PolNos);

  //elementos geom�tricos
  if(0) CriaElementos (4,4,*geomesh, IncidQuads);
  if(0) CriaElementos (4,8,*geomesh, Cubic3d);
  if(0) CriaElementos (1,4,*geomesh, Teste1Elems);
  if(0) CriaElementos (4,3,*geomesh, Teste2Elems);
  if(1) CriaElementos (1,4,*geomesh,  PolElems);

  //montagem de conectividades entre elementos
  geomesh->BuildConnectivity();

  //cria malha computacional
  TPZCompMesh *compmesh = new TPZCompMesh(geomesh);

  //cria material do problema
  TPZMaterial *placa = LerMaterial("placa.in",*compmesh);

  //transferindo para o material uma carga conhecida
  //placa->SetForcingFunction(LoadSolution);
  if(problema==1) (dynamic_cast<TPZPlaca *>(placa))->SetExactFunction(LoadSolution1);
  if(problema==2) (dynamic_cast<TPZPlaca *>(placa))->SetExactFunction(LoadSolution2);
  if(problema==3) (dynamic_cast<TPZPlaca *>(placa))->SetExactFunction(LoadSolution3);
  if(problema==11) (dynamic_cast<TPZPlaca *>(placa))->SetExactFunction(LoadSolution11);
  if(problema==12) (dynamic_cast<TPZPlaca *>(placa))->SetExactFunction(LoadSolution12);
  if(problema==13) (dynamic_cast<TPZPlaca *>(placa))->SetExactFunction(LoadSolution13);
  if(problema==14) (dynamic_cast<TPZPlaca *>(placa))->SetExactFunction(LoadSolution14);

  //cria condi��es de contorno
  if(0) CriaCondContTeste1(*geomesh);
  if(0) CriaCondContTeste2(*geomesh);
  if(1) CriaCondContTeste3(*geomesh);
  if(0) CriaCondContTeste4(*geomesh);

  //cria elementos computacionais
  compmesh->AutoBuild();

  if(1){
    //calculo da norma energia
    NormSol = SolutionError(compmesh,"Norma energia da solucao exata : ");
    //return 0;
  }

  //divide n�veis
  NivelDivide(compmesh);

  //ajusta elementos no contorno
  compmesh->AdjustBoundaryElements();

  //analysis do problema
  ExecutaAnalysis(*compmesh,placa);

  //calcula o erro energia da solu��o
  SolutionError(compmesh,"ERRO PERCENTUAL FINAL : ");

  //arquivo de saida de dados
  ofstream data("mesh.out");
  geomesh->Print(data);
  compmesh->Print(data);
  data.flush();
  data.close();
  erros.close();

  delete compmesh;
  delete geomesh;
  return 0;
}

/**
 * --> FIM PROGRAMA PRINCIPAL <-- --> FIM PROGRAMA PRINCIPAL <-- --> FIM PROGRAMA PRINCIPAL <--
 */

void CriaNos(int num, TPZGeoMesh &geomesh, double list [20][3] ){
  geomesh.NodeVec().Resize(num);
  TPZVec<REAL> coord(3);
  int i;
  if(problema==3){
    for(i=0;i<num;i++) {
      coord[0]= 2.0*list[i][0];
      coord[1]= 2.0*list[i][1];
      coord[2]= list[i][2];
      geomesh.NodeVec()[i].Initialize(coord,geomesh);
    }
  } else {
    for(i=0;i<num;i++) {
      coord[0]= list[i][0];
      coord[1]= list[i][1];
      coord[2]= list[i][2];
      geomesh.NodeVec()[i].Initialize(coord,geomesh);
    }
  }
}

void  CriaElementos(int numelem, int ncon, TPZGeoMesh &geomesh, int list[20][8]){

  TPZVec<int> index(ncon);
  int i,j,indice;
  for(i=0;i<numelem;i++) {
    for(j=0;j<ncon;j++) {
      index[j] = list [i][j];
    }
    if(ncon == 8) geomesh.CreateGeoElement(ECube,index,1,indice);
    if(ncon == 4) geomesh.CreateGeoElement(EQuadrilateral,index,1,indice);
    if(ncon == 3) geomesh.CreateGeoElement(ETriangle,index,1,indice);
  }
}

TPZMaterial *LerMaterial(char *filename, TPZCompMesh &cmesh) {

  ifstream input(filename);
  TPZFMatrix naxes(3,3);
  REAL ni1,ni2,h,E1,E2,G12,G13,G23,f;
  REAL n00,n01,n02,n10,n11,n12,n20,n21,n22;
  TPZVec<REAL> xf(6);
  int matindex;
  input >> matindex;
  input >> f   >>  h  >>
    E1  >> E2  >>
    G12 >> G13 >> G23 >>
    ni1 >> ni2;
  input >> n00 >> n01 >> n02;
  input >> n10 >> n11 >> n12;
  input >> n20 >> n21 >> n22;
  input >> xf[0] >> xf[1] >> xf[2] >> xf[3] >> xf[4] >> xf[5];
  naxes(0,0) =  n00;    naxes(0,1) =  n01;    naxes(0,2) =  n02;
  naxes(1,0) =  n10;    naxes(1,1) =  n11;    naxes(1,2) =  n12;
  naxes(2,0) =  n20;    naxes(2,1) =  n21;    naxes(2,2) =  n22;
  TPZMaterial *placa = new TPZPlaca(matindex,h,f,E1,E2,ni1,ni2,G12,G13,G23,naxes,xf);
  cmesh.InsertMaterialObject(placa);
  return placa;
}

void CriaCondContTeste3(TPZGeoMesh &gmesh){

  int indicematerial = 1;
  TPZAutoPointer<TPZMaterial> placa = gmesh.Reference()->FindMaterial(indicematerial);
  if(!placa){
    cout << "main::CriaCond material nao existe, CC nao criadas\n";
    cout << "\t\tindice material pedido : " << indicematerial << endl;
    return;
  }

  TPZGeoEl *elg0 = gmesh.FindElement(0);

  //malha computacional
  TPZCompMesh *cmesh = gmesh.Reference();

  //BIG number
  //   REAL big = 1.e12;
   //valor das CC
   TPZFMatrix val1(6,6,0.),val2(6,1,0.);

   TPZGeoElBC(elg0,4,-1,gmesh);
   TPZGeoElBC(elg0,5,-1,gmesh);
   TPZGeoElBC(elg0,6,-1,gmesh);
   TPZGeoElBC(elg0,7,-1,gmesh);
   TPZBndCond *bc = placa->CreateBC(placa,-1,0,val1,val2);
   bc->SetForcingFunction(BCSolution);
   cmesh->InsertMaterialObject(bc);

}

void CriaCondContTeste4(TPZGeoMesh &gmesh){

  int indicematerial = 1;
  TPZAutoPointer<TPZMaterial> placa = gmesh.Reference()->FindMaterial(indicematerial);
  if(!placa){
    cout << "main::CriaCond material nao existe, CC nao criadas\n";
    cout << "\t\tindice material pedido : " << indicematerial << endl;
    return;
  }

  TPZGeoEl *elg0 = gmesh.FindElement(0);
  TPZGeoEl *elg1 = gmesh.FindElement(1);
  TPZGeoEl *elg2 = gmesh.FindElement(2);
  TPZGeoEl *elg3 = gmesh.FindElement(3);

  //malha computacional
  TPZCompMesh *cmesh = gmesh.Reference();

  //BIG number
  //REAL big = 1.e12;
  //valor das CC
  TPZFMatrix val1(6,6,0.),val2(6,1,0.);

  TPZGeoElBC(elg0,3,-1,gmesh);
  TPZGeoElBC(elg1,3,-1,gmesh);
  TPZGeoElBC(elg2,3,-1,gmesh);
  TPZGeoElBC(elg3,3,-1,gmesh);
  TPZBndCond *bc = placa->CreateBC(placa,-1,0,val1,val2);
  bc->SetForcingFunction(BCSolution);
  cmesh->InsertMaterialObject(bc);

}

void ExecutaAnalysis(TPZCompMesh &cmesh,TPZMaterial *mat){

  //calcula o bloco da malha
   cmesh.InitializeBlock();

   //define o tipo de matriz de rigidez
   //matriz cheia n�o singular: requer um solver direto
   //TPZFStructMatrix strm(&cmesh);
   TPZSkylineStructMatrix strm(&cmesh);

   //cria objeto analysis
   TPZAnalysis an(&cmesh);

   //passa a matriz de rigidez para o objeto analysis
   an.SetStructuralMatrix(strm);

   //cria objeto solver
   TPZStepSolver solver;

   //define a decomposi��o LU
   solver.SetDirect(ELDLt);//ECholesky, ELDLt, ELU

   //associa o solver com o objeto analysis
   an.SetSolver(solver);

   //an�lise adaptativa
   if(1){
     REAL erro_adm;
     int numiter,resolut;
     cout << "main::ExecutaAnalysis error admissivel %\n";
     cin >> erro_adm;
     //tol = .05; cout << endl;//5%
     cout << "main::ExecutaAnalysis numero de iteracoes -> \n";
     cin >> numiter;
     //numiter = 20; cout << endl;
     cout << "main::ExecutaAnalysis resolucao  -> 0\n";
     //cin >> resolut;
     resolut = 0;
     erro_adm /= 100.0;
     AdaptativeProcedure(erro_adm,numiter,resolut,mat,an);
   }

   if(0){
     //monta a global baseada nas matrizes elementares e o bloco da malha
     cout << "main::ExecutaAnalysis comeca assemble\n";
     an.Assemble();

     //resolve o sistema com a matriz cheia e o solver direto com decomposi��o LU
     cout << "main::ExecutaAnalysis comeca solve\n";
     an.Solve();
   }

   //print do analysis
   ofstream dados("analysis.out");
   an.Print("FEM SOLUTION ",dados);

   if(0){
     //p�s processamento do main
     TPZGeoMesh *gmesh = an.Mesh()->Reference();
     ProcessamentoLocal(*gmesh,dados);

     //p�s processamento do DX
     PosProcessamento(an);
   }

   dados.flush();
   dados.close();
}

void PosProcessamento(TPZAnalysis &an){

  TPZVec<std::string> scalar(4);
  scalar[0] = "Deslocz";//var = 4
  scalar[1] = "Mn1";//var = 5
  scalar[2] = "Mn2";//var = 6
  scalar[3] = "Mn1n2";//var = 7
  TPZVec<std::string> vetorial(0);
  TPZAutoPointer<TPZMaterial> mat = an.Mesh()->FindMaterial(1);
  int dim = mat->Dimension();

  an.DefineGraphMesh(dim, scalar, vetorial, "PLACA.dx");
  cout << "\nmain::PosProcessamento arquivo de saida -> PLACA.dx\n";
  an.PostProcess(0);
  an.CloseGraphMesh();
}

void ProcessamentoLocal(TPZGeoMesh &gmesh,std::ostream &out) {

  int numpoints;
  int nel = gmesh.Reference()->ElementVec().NElements();
  int count = -1;
  while(count++<nel){
    for(int iel=0;iel<nel;iel++) {
      if(!gmesh.Reference()->ElementVec()[iel]) continue;
      int elemtype = gmesh.Reference()->ElementVec()[iel]->Type();
      if(elemtype==0 && elemtype==1) continue;
      TPZCompEl *el = gmesh.Reference()->ElementVec()[iel];
      TPZGeoEl *gel = el->Reference();
      if(el && gel) {
	if(gel->Id()==count){
	  out << "Elemento " << el->Reference()->Id() << endl;;
	  TPZManVector<REAL> sol(1);
	  TPZVec<REAL> csi(3,0.),x(3);
	  ifstream in("pontos.in");
	  in >> numpoints;
	  for(int p=0;p<numpoints;p++) {
	    in >> csi[0] >> csi[1] >>  csi[2];
	    gel->X(csi,x);
	    int var = 4;
	    gel->Reference()->Solution(csi,var,sol);
	    out << "solucao em x    = " << x[0] << ' ' << x[1] << ' ' << x[2] << endl;
	    if(sol.NElements()>0)
	      out << "               u = " << sol[0] << endl;
	  }
	  in.close();
	} else {
	  continue;
	}
      }
    }
  }
}

void NivelDivide(TPZCompMesh *cmesh){

  TPZVec<int> csub(0);
  int nivel;
  //cout << "\nmain::Divisao todos os elementos da malha serao divididos!\n";
  cout << "\nEntre nivel inicial da malha -> ";
  //cin >> nivel;
  nivel = 0; cout << endl;
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

REAL SolutionError(TPZCompMesh *cmesh,char *title){

  REAL Energy_error = 0.0,
    error1 = 0.0,
    error2 = 0.0;

  //erro da aproximacao FEM
  TPZManVector<REAL> flux(9);
  TPZManVector<REAL,3> errors(3,0.);
#warning Taken out an important call!!!
  std::cout << __PRETTY_FUNCTION__ << " Program this again\n";
  cmesh->EvaluateError(Solution,errors);
  Energy_error = errors[0];
  error1 = errors[1];
  error2 = errors[2];
  if(strcmp(title,"0")){
    if(NormSol < 0)
      cout << "\n" << title << Energy_error << "\n\n";
    else {
      cout << "\n" << title << (Energy_error/NormSol)*100.0 << "\n\n";
    }
  }
  return Energy_error;
}


REAL NormSolExact(TPZCompMesh *cmesh,char *title){
  //Dada a solu�p conhecida a sua norma pode ser calculada a qualquer momento
  //NormSol = NormSolExact(compmesh,"Norma energia da solucao exata : ");// <-- chamada
  TPZFMatrix *mat = dynamic_cast<TPZFMatrix *>(cmesh->Block().Matrix());
  TPZFMatrix copymat(*mat);
  copymat.Zero();
  cmesh->Block().SetMatrix(&copymat);
  REAL norm = SolutionError(cmesh,title);
  cmesh->Block().SetMatrix(mat);
  return norm;
}

void BCSolution(TPZVec<REAL> &x,TPZVec<REAL> &result){

  TPZFMatrix deriv(2,6);
  Solution(x,result,deriv);
}

void Solution(TPZVec<REAL> &x,TPZVec<REAL> &result,TPZFMatrix &deriv){

  TPZFMatrix eixos(3,3,0.);
  eixos(0,0) = 1.0;
  eixos(1,1) = 1.0;
  eixos(2,2) = 1.0;
  TPZFMatrix u(6,1);
  if(problema==1) LoadSolution1(eixos,x,u,deriv);
  if(problema==2) LoadSolution2(eixos,x,u,deriv);
  if(problema==3) LoadSolution3(eixos,x,u,deriv);
  if(problema==11) LoadSolution11(eixos,x,u,deriv);
  if(problema==12) LoadSolution12(eixos,x,u,deriv);
  if(problema==13) LoadSolution13(eixos,x,u,deriv);
  if(problema==14) LoadSolution14(eixos,x,u,deriv);
  for(int i=0;i<6;i++) result[i] = u(i,0);
}

void LoadSolution3(TPZFMatrix &axes,TPZVec<REAL> &X,TPZFMatrix &u,TPZFMatrix &du){

  REAL k,eps,div;
  static int key = 1;
  static REAL val;

  k = 0.0;//0.00001;
  eps = 0.4;//-1.0 � interessante
  div = 100.0;

  if(key){
    cout << "main::LoadSolution valor de div = ";
    cin >> val;
    cout << endl;
    key = 0;
  }
  div = val;

  REAL x = X[0],y=X[1];

  REAL x2 = (x-1.0)*(x-1.0);
  REAL y2 = (y-1.0)*(y-1.0);
  REAL eps2 = eps*eps;

  REAL expx = exp(-x2/eps2);
  REAL expy = exp(-y2/eps2);

  REAL exdivmk = expx/div+k;
  REAL eydivmk = expy/div+k;

  REAL exy = expx*expy;

  REAL exmeydivmk = 2.*expx*eydivmk / div / eps2;
  REAL eymexdivmk = 2.*expy*exdivmk / div / eps2;
  //u exact
  u(0,0) =  0.0;//u
  u(1,0) =  0.0;//v
  u(2,0) =  (expx/div+k)*(expy/div+k);//w
  u(3,0) = -( 2.*expy*(expx/div+k)*(y-1.) ) / div / eps2;//�x
  u(4,0) =  ( 2.*expx*(expy/div+k)*(x-1.) ) / div / eps2;//�y
  u(5,0) =  0.0;//�z

  //du exact
  TPZFMatrix dur(2,6);
  dur(0,0)  =  0.0;//du/dx
  dur(1,0)  =  0.0;//du/dy

  dur(0,1)  =  0.0;//dv/dx
  dur(1,1)  =  0.0;//dv/dx

  dur(0,2)  =  -exmeydivmk*(x-1.);//dw/dx
  dur(1,2)  =  -eymexdivmk*(y-1.);//dw/dy

  dur(0,3)  =  4.*exy*(x-1.)*(y-1.) / (div*div) / (eps2*eps2);//d�x/dx
  dur(1,3)  = -eymexdivmk + 2.*eymexdivmk*(y-1.)*(y-1.) / eps2;//d�x/dy

  dur(0,4)  =  exmeydivmk - 2.*exmeydivmk*(x-1.)*(x-1.) / eps2;//d�y/dx
  dur(1,4)  = -4.*exy*(x-1.)*(y-1.) / (div*div) / (eps2*eps2);//d�y/dy

  dur(0,5)  = 0.0;//d�z/dx
  dur(1,5)  = 0.0;//d�z/dy

  du(0,0) = axes(0,0)*dur(0,0) + axes(1,0)*dur(1,0);
  du(1,0) = axes(0,1)*dur(0,0) + axes(1,1)*dur(1,0);

  du(0,1) = axes(0,0)*dur(0,1) + axes(1,0)*dur(1,1);
  du(1,1) = axes(0,1)*dur(0,1) + axes(1,1)*dur(1,1);

  du(0,2) = axes(0,0)*dur(0,2) + axes(1,0)*dur(1,2);
  du(1,2) = axes(0,1)*dur(0,2) + axes(1,1)*dur(1,2);

  du(0,3) = axes(0,0)*dur(0,3) + axes(1,0)*dur(1,3);
  du(1,3) = axes(0,1)*dur(0,3) + axes(1,1)*dur(1,3);

  du(0,4) = axes(0,0)*dur(0,4) + axes(1,0)*dur(1,4);
  du(1,4) = axes(0,1)*dur(0,4) + axes(1,1)*dur(1,4);

  du(0,5) = axes(0,0)*dur(0,5) + axes(1,0)*dur(1,5);
  du(1,5) = axes(0,1)*dur(0,5) + axes(1,1)*dur(1,5);

}

void LoadSolution2(TPZFMatrix &axes,TPZVec<REAL> &X,TPZFMatrix &u,TPZFMatrix &du){

  REAL x = X[0];
  REAL y = X[1];

  REAL c = 1.0;
  REAL p=1.0,q=1.0;

  REAL xp = pow(x,p);
  REAL yq = pow(y,q);

  REAL xpm1 = pow(x,p-1.);
  REAL yqm1 = pow(y,q-1.);

  REAL xm1p = pow(x-1.,p);
  REAL ym1q = pow(y-1.,q);

  REAL xm1pm1 = pow(x-1.,p-1.);
  REAL ym1qm1 = pow(y-1.,q-1.);

  REAL yqm2 = pow(y,q-2.);
  REAL ym1qm2 = pow(y-1.,q-2.);
  REAL xpm2 = pow(x,p-2.);
  REAL xm1pm2 = pow(x-1,p-2.);


  u(0,0) =  0.;
  u(1,0) =  0.;
  u(2,0) = -c*xm1p*xp*ym1q*yq;
  u(3,0) = -c*q*xm1p*xp*ym1q*yqm1 + c*q*xm1p*xp*ym1qm1*yq;
  u(4,0) =  c*p*xm1p*xpm1*ym1q*yq - c*p*xm1pm1*xp*ym1q*yq;
  u(5,0) =  0.;

  //du exact
  TPZFMatrix dur(2,6);
  dur(0,0) =  0.;//dudx
  dur(1,0) =  0.;//dudy
  dur(0,1) =  0.;
  dur(1,1) =  0.;
  dur(0,2) = -c*p*xm1p*xpm1*ym1q*yq + c*p*xm1pm1*xp*ym1q*yq;
  dur(1,2) = -c*q*xm1p*xp*ym1q*yqm1 + c*q*xm1p*xp*ym1qm1*yq;
  dur(0,3) = -c*p*q*xm1p*xpm1*ym1q*yqm1 + c*p*q*xm1pm1*xp*ym1q*yqm1 + c*p*q*xm1p*xpm1*ym1qm1*yq + c*p*q*xm1pm1*xp*ym1qm1*yq;
  dur(1,3) = -c*(q-1)*q*xm1p*xp*ym1q*yqm2 + 2*c*q*q*xm1p*xp*ym1qm1*yqm1 + c*(q-1)*q*xm1p*xp*ym1qm2*yq;
  dur(0,4) =  c*(p-1)*p*xm1p*xpm2*ym1q*yq - 2*c*p*p*xm1pm1*xpm1*ym1q*yq - c*(p-1)*p*xm1pm2*xp*ym1q*yq;
  dur(1,4) =  c*p*q*xm1p*xpm1*ym1q*yqm1 - c*p*q*xm1pm1*xp*ym1q*yqm1 - c*p*q*xm1p*xpm1*ym1qm1*yq - c*p*q*xm1pm1*xp*ym1qm1*yq;
  dur(0,5) =  0.;
  dur(1,5) =  0.;

  du(0,0) = axes(0,0)*dur(0,0) + axes(1,0)*dur(1,0);
  du(1,0) = axes(0,1)*dur(0,0) + axes(1,1)*dur(1,0);

  du(0,1) = axes(0,0)*dur(0,1) + axes(1,0)*dur(1,1);
  du(1,1) = axes(0,1)*dur(0,1) + axes(1,1)*dur(1,1);

  du(0,2) = axes(0,0)*dur(0,2) + axes(1,0)*dur(1,2);
  du(1,2) = axes(0,1)*dur(0,2) + axes(1,1)*dur(1,2);

  du(0,3) = axes(0,0)*dur(0,3) + axes(1,0)*dur(1,3);
  du(1,3) = axes(0,1)*dur(0,3) + axes(1,1)*dur(1,3);

  du(0,4) = axes(0,0)*dur(0,4) + axes(1,0)*dur(1,4);
  du(1,4) = axes(0,1)*dur(0,4) + axes(1,1)*dur(1,4);

  du(0,5) = axes(0,0)*dur(0,5) + axes(1,0)*dur(1,5);
  du(1,5) = axes(0,1)*dur(0,5) + axes(1,1)*dur(1,5);

}

void LoadSolution11(TPZFMatrix &axes,TPZVec<REAL> &x,TPZFMatrix &u,TPZFMatrix &du){
  //u exact
  u(0,0) =  0.0;
  u(1,0) =  0.0;
  u(2,0) = -x[0]*x[0]*x[1]*x[1];//x*y
  u(3,0) = -2.0*x[0]*x[0]*x[1];//dw/dy
  u(4,0) =  2.0*x[0]*x[1]*x[1];//-dw/dx
  u(5,0) =  0.0;//0

  //du exact
  TPZFMatrix dur(2,6);
  dur(0,0) =  0.0;
  dur(1,0) =  0.0;
  dur(0,1) =  0.0;
  dur(1,1) =  0.0;
  dur(0,2) = -2.0*x[0]*x[1]*x[1];
  dur(1,2) = -2.0*x[0]*x[0]*x[1];
  dur(0,3) = -4.0*x[0]*x[1];
  dur(1,3) = -2.0*x[0]*x[0];
  dur(0,4) =  2.0*x[1]*x[1];
  dur(1,4) =  4.0*x[0]*x[1];
  dur(0,5) =  0.0;
  dur(1,5) =  0.0;

  du(0,0) = axes(0,0)*dur(0,0) + axes(1,0)*dur(1,0);
  du(1,0) = axes(0,1)*dur(0,0) + axes(1,1)*dur(1,0);

  du(0,1) = axes(0,0)*dur(0,1) + axes(1,0)*dur(1,1);
  du(1,1) = axes(0,1)*dur(0,1) + axes(1,1)*dur(1,1);

  du(0,2) = axes(0,0)*dur(0,2) + axes(1,0)*dur(1,2);
  du(1,2) = axes(0,1)*dur(0,2) + axes(1,1)*dur(1,2);

  du(0,3) = axes(0,0)*dur(0,3) + axes(1,0)*dur(1,3);
  du(1,3) = axes(0,1)*dur(0,3) + axes(1,1)*dur(1,3);

  du(0,4) = axes(0,0)*dur(0,4) + axes(1,0)*dur(1,4);
  du(1,4) = axes(0,1)*dur(0,4) + axes(1,1)*dur(1,4);

  du(0,5) = axes(0,0)*dur(0,5) + axes(1,0)*dur(1,5);
  du(1,5) = axes(0,1)*dur(0,5) + axes(1,1)*dur(1,5);
}

void LoadSolution10(TPZFMatrix &axes,TPZVec<REAL> &x,TPZFMatrix &u,TPZFMatrix &du){

  //u exact
  u(0,0) =  0.0;
  u(1,0) =  0.0;
  u(2,0) = -x[0]*x[0]*x[1]*x[1]*x[1];//x*y
  u(3,0) = -3.0*x[0]*x[0]*x[1]*x[1];//dw/dy
  u(4,0) =  2.0*x[0]*x[1]*x[1]*x[1];//-dw/dx
  u(5,0) =  0.0;//0

  //du exact
  TPZFMatrix dur(2,6);
  dur(0,0) =  0.0;
  dur(1,0) =  0.0;
  dur(0,1) =  0.0;
  dur(1,1) =  0.0;
  dur(0,2) = -2.0*x[0]*x[1]*x[1]*x[1];
  dur(1,2) = -3.0*x[0]*x[0]*x[1]*x[1];
  dur(0,3) = -6.0*x[0]*x[1]*x[1];
  dur(1,3) = -6.0*x[0]*x[0]*x[1];
  dur(0,4) =  2.0*x[1]*x[1]*x[1];
  dur(1,4) =  6.0*x[0]*x[1]*x[1];
  dur(0,5) =  0.0;
  dur(1,5) =  0.0;

  du(0,0) = axes(0,0)*dur(0,0) + axes(1,0)*dur(1,0);
  du(1,0) = axes(0,1)*dur(0,0) + axes(1,1)*dur(1,0);

  du(0,1) = axes(0,0)*dur(0,1) + axes(1,0)*dur(1,1);
  du(1,1) = axes(0,1)*dur(0,1) + axes(1,1)*dur(1,1);

  du(0,2) = axes(0,0)*dur(0,2) + axes(1,0)*dur(1,2);
  du(1,2) = axes(0,1)*dur(0,2) + axes(1,1)*dur(1,2);

  du(0,3) = axes(0,0)*dur(0,3) + axes(1,0)*dur(1,3);
  du(1,3) = axes(0,1)*dur(0,3) + axes(1,1)*dur(1,3);

  du(0,4) = axes(0,0)*dur(0,4) + axes(1,0)*dur(1,4);
  du(1,4) = axes(0,1)*dur(0,4) + axes(1,1)*dur(1,4);

  du(0,5) = axes(0,0)*dur(0,5) + axes(1,0)*dur(1,5);
  du(1,5) = axes(0,1)*dur(0,5) + axes(1,1)*dur(1,5);
}

void LoadSolution1(TPZFMatrix &axes,TPZVec<REAL> &x,TPZFMatrix &u,TPZFMatrix &du){

  //u exact
  u(0,0) =  0.0;
  u(1,0) =  0.0;
  u(2,0) = -x[0]*x[1];//x*y
  u(3,0) = -x[0];//dw/dy
  u(4,0) =  x[1];//-dw/dx
  u(5,0) =  0.0;//0

  //du exact
  TPZFMatrix dur(2,6);
  dur(0,0) =  0.0;
  dur(1,0) =  0.0;
  dur(0,1) =  0.0;
  dur(1,1) =  0.0;
  dur(0,2) = -x[1];
  dur(1,2) = -x[0];
  dur(0,3) = -1.0;
  dur(1,3) =  0.0;
  dur(0,4) =  0.0;
  dur(1,4) =  1.0;
  dur(0,5) =  0.0;
  dur(1,5) =  0.0;

  du(0,0) = axes(0,0)*dur(0,0) + axes(1,0)*dur(1,0);
  du(1,0) = axes(0,1)*dur(0,0) + axes(1,1)*dur(1,0);

  du(0,1) = axes(0,0)*dur(0,1) + axes(1,0)*dur(1,1);
  du(1,1) = axes(0,1)*dur(0,1) + axes(1,1)*dur(1,1);

  du(0,2) = axes(0,0)*dur(0,2) + axes(1,0)*dur(1,2);
  du(1,2) = axes(0,1)*dur(0,2) + axes(1,1)*dur(1,2);

  du(0,3) = axes(0,0)*dur(0,3) + axes(1,0)*dur(1,3);
  du(1,3) = axes(0,1)*dur(0,3) + axes(1,1)*dur(1,3);

  du(0,4) = axes(0,0)*dur(0,4) + axes(1,0)*dur(1,4);
  du(1,4) = axes(0,1)*dur(0,4) + axes(1,1)*dur(1,4);

  du(0,5) = axes(0,0)*dur(0,5) + axes(1,0)*dur(1,5);
  du(1,5) = axes(0,1)*dur(0,5) + axes(1,1)*dur(1,5);
}

static clock_t fBegin=0,fInit=0;
void CoutTime(clock_t &start,const char *title);
void AdaptativeProcedure(REAL erro_adm,int numiter,int resolution,TPZMaterial *mat,TPZAnalysis &an) {
  cout << "PZAnalysis::IterativeProcessTest beginning of the iterative process, general time 0\n";
  int iter = 0;
  TPZCompMesh *cmesh = an.Mesh();
  int numeq = cmesh->NEquations();
  cout << "\nmain::AdaptativeProcedure numero de equacoes    : " << numeq;
  int numel = NumElComp(*cmesh);
  fBegin = clock();
  an.Solution().Zero();
  an.Run();
  an.LoadSolution();
  cout << "TPZIterativeAnalysis::AdaptativeProcedure iteracao = " << ++iter << endl;
  CoutTime(fBegin,"TPZIterativeAnalysis::AdaptativeProcedure end system solution first iteration");
  fBegin = clock();
  SolutionError(cmesh,"Erro Energia  %  : ");
  CoutTime(fBegin,"TPZIterativeAnalysis::AdaptativeProcedure end solution error");
  REAL normsol = Norm(an.Solution());
  REAL tolmax = 1.e15;
  REAL tol =  erro_adm * NormSol / sqrt(REAL(numel));
  fBegin = clock();
  int refined = AdaptiveMesh(*cmesh,tol,iter);
  CoutTime(fBegin,"TPZIterativeAnalysis::AdaptativeProcedure end adaptative mesh");

  while( iter < numiter && normsol < tolmax && refined > 0 ) {

    numeq = cmesh->NEquations();
    cout << "\nmain::AdaptativeProcedure numero de equacoes    : " << numeq;
    numel = NumElComp(*cmesh);
    fBegin = clock();
    an.SetBlockNumber();
    an.Solution().Zero();
    an.Run();
    an.LoadSolution();
    normsol = Norm(an.Solution());
    CoutTime(fBegin,"TPZIterativeAnalysis::AdaptativeProcedure end system solution actual iteration");
    CoutTime(fInit,"TPZIterativeAnalysis::AdaptativeProcedure accumulated time");
    fBegin = clock();
    SolutionError(cmesh,"Erro Energia % : ");
    CoutTime(fBegin,"TPZIterativeAnalysis::AdaptativeProcedure end solution error");
    tol =  erro_adm * NormSol / sqrt(REAL(numel));
    fBegin = clock();
    refined = AdaptiveMesh(*cmesh,tol,iter);
    CoutTime(fBegin,"TPZIterativeAnalysis::AdaptativeProcedure end adaptative mesh");
    cout << "TPZIterativeAnalysis::AdaptativeProcedure iteracao = " << ++iter << endl;

  }

  if( normsol > tolmax ){
    cout << "\nTPZIterativeAnalysis::AdaptativeProcedure the iterative process stopped due the great norm "
	 << "of the solution, norm solution = " << normsol << endl;
  }
  CoutTime(fInit,"TPZIterativeAnalysis::AdaptativeProcedure general time of iterative process");
  if(!refined && iter <= numiter){
    cout << "\nTPZIterativeAnalysis::AdaptativeProcedure \n\n"
	 <<"\t\t\t\t\t<-> TOLERANCIA ATINGIDA <->\n\n";
  }
  if(refined == -1) cout << "\n\nPZIterativeAnalysis::AdaptativeProcedure STOP FOR MAXIM DOF\n\n";
  TPZVec<std::string> scalar(2),vector(0);
  scalar[0] = "Deslocz";
  scalar[1] = "POrder";
  //scalar[1] = "Mn1";
  //scalar[2] = "Mn2";
  int dim = mat->Dimension();

  an.DefineGraphMesh(dim, scalar, vector, "PLACA.dx");
  cout << "\nmain::PosProcessamento arquivo de saida -> PLACA.dx\n";
  an.SetTime(0.1);
  an.PostProcess(resolution);
  cout << "\nTPZIterativeAnalysis:: out file : PLACA.dx\n";
  an.CloseGraphMesh();
}

void CoutTime(clock_t &start,const char *title){
    clock_t end = clock();
    cout << title <<  endl;
    clock_t segundos = ((end - start)/CLOCKS_PER_SEC);
    cout << segundos << " segundos" << endl;
    cout << segundos/60.0 << " minutos" << endl << endl;
}

int AdaptiveMesh(TPZCompMesh &cmesh,REAL tol,int iter){

  int numdivide = 0;
  erros << "\n\n\n *-**-**-*\tmain::AdaptiveMesh ITERACAO = " << iter;
  erros << "\t*-**-**-*\n\n";
  int nel = cmesh.ElementVec().NElements();
  int numeq = cmesh.NEquations();
  if(numeq > 5000) return -1;
  int i,ordmax = 5;
  TPZBlock *flux;
  REAL energy,error2,error3;
  TPZVec<int> sub(0);
  TPZStack<int> nindex,pindex;
  nindex.Resize(0);
  pindex.Resize(0);
  TPZCompEl *com;
  TPZInterpolatedElement *intel;
  REAL maxerror = 0.0,minerror=1000.0;
  int elmaxerr,elid,order;
  REAL MeshError = SolutionError(&cmesh,"0");
  REAL user_error = tol*sqrt(REAL(NumElComp(cmesh))) / NormSol;
  for(i=0;i<nel;i++){
    com = cmesh.ElementVec()[i];
    if(!com) continue;
    //if(com->Type() == 16) continue
    if(com->Material()->Id() < 0) continue;
    intel = dynamic_cast<TPZInterpolatedElement *>(com);
    TPZManVector<REAL,3> errors(3,0.);
    intel->EvaluateError(Solution,errors,flux);
    energy = errors[0];
    error2 = errors[1];
    error3 = errors[2];
//    intel->EvaluateError(Solution,energy,error2,flux,error3);
    elid = com->Reference()->Id();
    erros << "\nElemento id   = " << elid;
    erros << "\nNorma energia absoluta = " << energy;
    if(maxerror < energy){
      maxerror = energy;
      elmaxerr = elid;
    }
    if(minerror > energy){
      minerror = energy;
    }
    intel = dynamic_cast<TPZInterpolatedElement *>(com);
    order = intel->SideOrder(intel->NConnects()-1);
    if(energy > tol && nel < 65){//251, 601
      nindex.Push(i);
      numdivide++;
    }
    else if(energy > tol && order < ordmax){
      pindex.Push(i);
      numdivide++;
    }
  }
  erros.flush();
  if(maxerror){
    cout << "main::AdaptiveMesh Maximo erro absoluto da malha : " << maxerror;
    cout << "\nmain::AdaptiveMesh elemento de id = " << elid << endl;
  }
  //erro da malha e maximo erro dos elementos s�o menor que a tolerancia para o erro global
  if(maxerror < user_error && MeshError < user_error) return 0;
  int nsize = nindex.NElements();
  int psize = pindex.NElements();//,neworder;
  if( maxerror > tol && !numdivide ) return -1;//erro uniforme n�o atingido, n�o h� adapta��o
  else if( maxerror <= tol ) return 0;//erro uniforme atingido
  for(i=0;i<nsize;i++){
    com = cmesh.ElementVec()[nindex[i]];
    cmesh.Divide(nindex[i],sub,1);
    numdivide++;
  }
  for(i=0;i<psize;i++){
     com = cmesh.ElementVec()[pindex[i]];
     intel = dynamic_cast<TPZInterpolatedElement *>(com);
     order = intel->SideOrder(intel->NConnects()-1);
//     neworder = order < ordmax ? (order +1) : order;
    //neworder = order <  Nivel(intel->Reference()) ? (order +1) : order;
//     if(neworder > order)
      intel->PRefine(order+1);
  }
  cmesh.AdjustBoundaryElements();
  return 1;
  //cmesh.InitializeBlock();
}

int NumElComp(TPZCompMesh &cmesh){

  int nel = cmesh.ElementVec().NElements();
  int i,num = 0,numbc = 0;
  for(i=0;i<nel;i++){
    TPZCompEl *comp = cmesh.ElementVec()[i];
    if(!comp) continue;
    if(comp->Material()->Id() < 0){
      numbc++;
      continue;
    }
    num++;
  }
  cout << "\nmain::NumElComp numero de elementos de volume   : " << num;
  cout << "\nmain::NumElComp numero de elementos de contorno : " << numbc << endl << endl;
  return num;
}


void LoadSolution13(TPZFMatrix &axes,TPZVec<REAL> &X,TPZFMatrix &u,TPZFMatrix &du){

  REAL x=X[0], y=X[1];
  //u exact
  u(0,0) =  0.0;
  u(1,0) =  0.0;
  u(2,0) =  y*y;//x2*y3
  u(3,0) =  2*y;//dw/dy
  u(4,0) =  0.0;//-dw/dx
  u(5,0) =  0.0;//0

  //du exact
  TPZFMatrix dur(2,6);
  dur(0,0) =  0.0;
  dur(1,0) =  0.0;
  dur(0,1) =  0.0;
  dur(1,1) =  0.0;
  dur(0,2) =  0.0;
  dur(1,2) =  2.0*y;
  dur(0,3) =  0.0;
  dur(1,3) =  2.0;
  dur(0,4) =  0.0;
  dur(1,4) =  0.0;
  dur(0,5) =  0.0;
  dur(1,5) =  0.0;

  du(0,0) = axes(0,0)*dur(0,0) + axes(1,0)*dur(1,0);
  du(1,0) = axes(0,1)*dur(0,0) + axes(1,1)*dur(1,0);

  du(0,1) = axes(0,0)*dur(0,1) + axes(1,0)*dur(1,1);
  du(1,1) = axes(0,1)*dur(0,1) + axes(1,1)*dur(1,1);

  du(0,2) = axes(0,0)*dur(0,2) + axes(1,0)*dur(1,2);
  du(1,2) = axes(0,1)*dur(0,2) + axes(1,1)*dur(1,2);

  du(0,3) = axes(0,0)*dur(0,3) + axes(1,0)*dur(1,3);
  du(1,3) = axes(0,1)*dur(0,3) + axes(1,1)*dur(1,3);

  du(0,4) = axes(0,0)*dur(0,4) + axes(1,0)*dur(1,4);
  du(1,4) = axes(0,1)*dur(0,4) + axes(1,1)*dur(1,4);

  du(0,5) = axes(0,0)*dur(0,5) + axes(1,0)*dur(1,5);
  du(1,5) = axes(0,1)*dur(0,5) + axes(1,1)*dur(1,5);
}



void LoadSolution14(TPZFMatrix &axes,TPZVec<REAL> &X,TPZFMatrix &u,TPZFMatrix &du){

  REAL x=X[0], y=X[1];
  //u exact
  u(0,0) =  0.0;
  u(1,0) =  0.0;
  u(2,0) =  x*x*y;//x2*y3
  u(3,0) =  x*x;//dw/dy
  u(4,0) =  -2.0*x*y;//-dw/dx
  u(5,0) =  0.0;//0

  //du exact
  TPZFMatrix dur(2,6);
  dur(0,0) =  0.0;
  dur(1,0) =  0.0;
  dur(0,1) =  0.0;
  dur(1,1) =  0.0;
  dur(0,2) =  2.0*x*y;
  dur(1,2) =  x*x;
  dur(0,3) =  2.0*x;
  dur(1,3) =  0.0;
  dur(0,4) = -2.0*y;
  dur(1,4) = -2.0*x;
  dur(0,5) =  0.0;
  dur(1,5) =  0.0;

  du(0,0) = axes(0,0)*dur(0,0) + axes(1,0)*dur(1,0);
  du(1,0) = axes(0,1)*dur(0,0) + axes(1,1)*dur(1,0);

  du(0,1) = axes(0,0)*dur(0,1) + axes(1,0)*dur(1,1);
  du(1,1) = axes(0,1)*dur(0,1) + axes(1,1)*dur(1,1);

  du(0,2) = axes(0,0)*dur(0,2) + axes(1,0)*dur(1,2);
  du(1,2) = axes(0,1)*dur(0,2) + axes(1,1)*dur(1,2);

  du(0,3) = axes(0,0)*dur(0,3) + axes(1,0)*dur(1,3);
  du(1,3) = axes(0,1)*dur(0,3) + axes(1,1)*dur(1,3);

  du(0,4) = axes(0,0)*dur(0,4) + axes(1,0)*dur(1,4);
  du(1,4) = axes(0,1)*dur(0,4) + axes(1,1)*dur(1,4);

  du(0,5) = axes(0,0)*dur(0,5) + axes(1,0)*dur(1,5);
  du(1,5) = axes(0,1)*dur(0,5) + axes(1,1)*dur(1,5);
}

void LoadSolution12(TPZFMatrix &axes,TPZVec<REAL> &X,TPZFMatrix &u,TPZFMatrix &du){

  REAL x=X[0], y=X[1];
  //u exact
  u(0,0) =  0.0;
  u(1,0) =  0.0;
  u(2,0) =  x*x;//x2*y3
  u(3,0) =  0;//dw/dy
  u(4,0) =  -2.0*x;//-dw/dx
  u(5,0) =  0.0;//0

  //du exact
  TPZFMatrix dur(2,6);
  dur(0,0) =  0.0;
  dur(1,0) =  0.0;
  dur(0,1) =  0.0;
  dur(1,1) =  0.0;
  dur(0,2) =  2.0*x;
  dur(1,2) =  0.0;
  dur(0,3) =  0.0;
  dur(1,3) =  0.0;
  dur(0,4) = -2.0;
  dur(1,4) =  0.0;
  dur(0,5) =  0.0;
  dur(1,5) =  0.0;

  du(0,0) = axes(0,0)*dur(0,0) + axes(1,0)*dur(1,0);
  du(1,0) = axes(0,1)*dur(0,0) + axes(1,1)*dur(1,0);

  du(0,1) = axes(0,0)*dur(0,1) + axes(1,0)*dur(1,1);
  du(1,1) = axes(0,1)*dur(0,1) + axes(1,1)*dur(1,1);

  du(0,2) = axes(0,0)*dur(0,2) + axes(1,0)*dur(1,2);
  du(1,2) = axes(0,1)*dur(0,2) + axes(1,1)*dur(1,2);

  du(0,3) = axes(0,0)*dur(0,3) + axes(1,0)*dur(1,3);
  du(1,3) = axes(0,1)*dur(0,3) + axes(1,1)*dur(1,3);

  du(0,4) = axes(0,0)*dur(0,4) + axes(1,0)*dur(1,4);
  du(1,4) = axes(0,1)*dur(0,4) + axes(1,1)*dur(1,4);

  du(0,5) = axes(0,0)*dur(0,5) + axes(1,0)*dur(1,5);
  du(1,5) = axes(0,1)*dur(0,5) + axes(1,1)*dur(1,5);
}

