
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
#include "pzelgc3d.h"
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

#include <math.h>
using namespace std;

//#include "BCTest.cc"

//static double quadrilatero[4][3]  = { {0.,0.,0.},{1.,0.,0.},{0.,1.,0.},{1.,1.,0.} };

static double NosQuads[9][3]  = { {-1.,-1.,0.},{0.,-1.,0.},{1.,-1.,0.},{-1.,0.,0.},{0.,0.,0.},
				  {1.,0.,0.},{-1.,1.,0.},{0.,1.,0.},{1.,1.,0.} };
static int    IncidQuads[4][8] = {{0,1,4,3},{1,2,5,4},{4,5,8,7},{3,4,7,6}};

static double Teste1Nos[4][3]  = { {0.,0.,0.},{2.,0.,0.},{2.,2.,0.},{0.,2.,0.}};
static int    Teste1Elems[1][8] = {{0,1,2,3}};


static double Teste2Nos[5][3]  = { {0.,0.,0.},{2.,0.,0.},{2.,2.,0.},{0.,2.,0.},{1.,1.,0.}};
static int    Teste2Elems[4][8] = {{0,1,4},{1,2,4},{2,3,4},{3,0,4}};


static double Nodos3d[18][3] = {{-1.,-1.,0.},{0.,-1.,0.},{1.,-1.,0.},{-1.,0.,0.,},{0.,0.,0.},
				{1.,0.,0.},{-1.,1.,0.},{0.,1.,0.},{1.,1.,0.},
				{-1.,-1.,.2},{0.,-1.,.2},{1.,-1.,.2},{-1.,0.,.2},{0.,0.,.2},
				{1.,0.,.2},{-1.,1.,.2},{0.,1.,.2},{1.,1.,.2}};

static int  Cubic3d[4][8]  = { {0,1,4,3,9,10,13,12},{1,2,5,4,10,11,14,13},
			       {4,5,8,7,13,14,17,16},{3,4,7,6,12,13,16,15} };

void ProcessamentoLocal(TPZGeoMesh &gmesh,ostream &out);
void CriaNos(int num, TPZGeoMesh &geomesh, double list [20][3]);
void CriaElementos(int num,int size,TPZGeoMesh &geomesh, int list[20][8] );
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
void BCSolution(TPZVec<REAL> &x,TPZVec<REAL> &result);
void Solution(TPZVec<REAL> &x,TPZVec<REAL> &result,TPZFMatrix &deriv);
REAL SolutionError(TPZCompMesh *cmesh);
TPZMaterial *LerMaterial(char *filename,TPZCompMesh &cmesh);
void AdaptativeProcedure(REAL error,int numiter,int marcha,int resolution,TPZMaterial *mat,TPZAnalysis &an);
void AdaptativeProcedure2(REAL tol,int numiter,int marcha,int resolution,TPZMaterial *mat,TPZAnalysis &an);
int AdaptiveMesh(TPZCompMesh &cmesh,REAL tol,int iter);
int NumElComp(TPZCompMesh &cmesh);
ofstream erros("ERROS.out");
/**
 * -> PROGRAMA PRINCIPAL <- -> PROGRAMA PRINCIPAL <- -> PROGRAMA PRINCIPAL <- -> PROGRAMA PRINCIPAL <-
 */

int main() {



  //inserindo a ordem de interpola� dos elementos e do espao
  int ord;
  cout << "Entre ordem 1,2,3,4,5 : -> 1";
  //cin >> ord;
  ord = 1; cout << endl;
  TPZCompEl::SetgOrder(ord);

  //malha geometrica
  TPZGeoMesh *geomesh = new TPZGeoMesh;

  geomesh->Reference();

  //cria n� geom�ricos
  if(0) CriaNos(9, *geomesh, NosQuads);
  if(0) CriaNos(18, *geomesh, Nodos3d);
  if(0) CriaNos(4, *geomesh, Teste1Nos);
  if(1) CriaNos(5, *geomesh, Teste2Nos);

  //elementos geom�ricos
  if(0) CriaElementos (4,4,*geomesh, IncidQuads);
  if(0) CriaElementos (4,8,*geomesh, Cubic3d);
  if(0) CriaElementos (1,4,*geomesh, Teste1Elems);
  if(1) CriaElementos (4,3,*geomesh, Teste2Elems);

  //montagem de conectividades entre elementos
  geomesh->BuildConnectivity2();

  //cria malha computacional
  TPZCompMesh *compmesh = new TPZCompMesh(geomesh);

  //cria material do problema
  TPZMaterial *placa = LerMaterial("placa.in",*compmesh);

  //transferindo para o material uma carga conhecida
  //placa->SetForcingFunction(LoadSolution);
  (dynamic_cast<TPZPlaca *>(placa))->SetExactFunction(LoadSolution);

  //cria condi�s de contorno
  if(0) CriaCondContTeste1(*geomesh);
  if(0) CriaCondContTeste2(*geomesh);
  if(0) CriaCondContTeste3(*geomesh);
  if(1) CriaCondContTeste4(*geomesh);

  //cria elementos computacionais
  compmesh->AutoBuild();

  //divide n�eis
  NivelDivide(compmesh);

  //ajusta elementos no contorno
  compmesh->AdjustBoundaryElements();

  //analysis do problema
  ExecutaAnalysis(*compmesh,placa);

  //calcula o erro energia da solu�
  SolutionError(compmesh);

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
  for(i=0;i<num;i++) {
    coord[0]= list[i][0];
    coord[1]= list[i][1];
    coord[2]= list[i][2];
    geomesh.NodeVec()[i].Initialize(coord,geomesh);

  }
}

void  CriaElementos(int num, int size, TPZGeoMesh &geomesh, int list[20][8]){

  TPZVec<int> index(size);
  int i,j,indice;
  for(i=0;i<num;i++) {
    for(j=0;j<size;j++) {
      index[j] = list [i][j];
    }
    if(size == 8) geomesh.CreateGeoElement(ECube,index,1,indice);
    if(size == 4) geomesh.CreateGeoElement(EQuadrilateral,index,1,indice);
    if(size == 3) geomesh.CreateGeoElement(ETriangle,index,1,indice);
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
  TPZMaterial *placa = gmesh.Reference()->FindMaterial(indicematerial);
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
   TPZBndCond *bc = placa->CreateBC(-1,0,val1,val2);
   bc->SetForcingFunction(BCSolution);
   cmesh->InsertMaterialObject(bc);

}

void CriaCondContTeste4(TPZGeoMesh &gmesh){

  int indicematerial = 1;
  TPZMaterial *placa = gmesh.Reference()->FindMaterial(indicematerial);
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
  TPZBndCond *bc = placa->CreateBC(-1,0,val1,val2);
  bc->SetForcingFunction(BCSolution);
  cmesh->InsertMaterialObject(bc);

}

void ExecutaAnalysis(TPZCompMesh &cmesh,TPZMaterial *mat){

  //calcula o bloco da malha
   cmesh.InitializeBlock();

   //define o tipo de matriz de rigidez
   //matriz cheia n� singular: requer um solver direto
   //TPZFStructMatrix strm(&cmesh);
   TPZSkylineStructMatrix strm(&cmesh);

   //cria objeto analysis
   TPZAnalysis an(&cmesh);

   //passa a matriz de rigidez para o objeto analysis
   an.SetStructuralMatrix(strm);

   //cria objeto solver
   TPZStepSolver solver;

   //define a decomposi� LU
   solver.SetDirect(ELDLt);//ECholesky, ELDLt, ELU

   //associa o solver com o objeto analysis
   an.SetSolver(solver);

   //an�ise adaptativa
   if(1){
     REAL tol;
     int numiter,marcha,resolut;
     cout << "main::ExecutaAnalysis tolerancia \n";
     cin >> tol;
     //tol = .05; cout << endl;//5%
     cout << "main::ExecutaAnalysis numero de iteracoes -> \n";
     //cin >> numiter;
     numiter = 20; cout << endl;
     cout << "main::ExecutaAnalysis marcha     -> \n";
     //cin >> marcha;
     marcha = 1;
     cout << "main::ExecutaAnalysis resolucao  -> 0\n";
     //cin >> resolut;
     resolut = 0;

     AdaptativeProcedure(tol,numiter,marcha,resolut,mat,an);
   }

   if(0){
     //monta a global baseada nas matrizes elementares e o bloco da malha
     cout << "main::ExecutaAnalysis comeca assemble\n";
     an.Assemble();

     //resolve o sistema com a matriz cheia e o solver direto com decomposi� LU
     cout << "main::ExecutaAnalysis comeca solve\n";
     an.Solve();
   }

   //print do analysis
   ofstream dados("analysis.out");
   an.Print("FEM SOLUTION ",dados);

   if(0){
     //p� processamento do main
     TPZGeoMesh *gmesh = an.Mesh()->Reference();
     ProcessamentoLocal(*gmesh,dados);

     //p� processamento do DX
     PosProcessamento(an);
   }

   dados.flush();
   dados.close();
}

void PosProcessamento(TPZAnalysis &an){

  TPZVec<char *> scalar(4);
  scalar[0] = "Deslocz";//var = 4
  scalar[1] = "Mn1";//var = 5
  scalar[2] = "Mn2";//var = 6
  scalar[3] = "Mn1n2";//var = 7
  TPZVec<char *> vetorial(0);
  TPZMaterial *mat = an.Mesh()->FindMaterial(1);
  int dim = mat->Dimension();
  TPZDXGraphMesh graph(an.Mesh(),dim,mat,scalar,vetorial);
  ofstream *dxout = new ofstream("PLACA.dx"); //visualiza com o DX
  cout << "\nmain::PosProcessamento arquivo de saida -> PLACA.dx\n";
  graph.SetOutFile(*dxout);
  graph.SetResolution(0);
  graph.DrawMesh(dim);
  REAL time = 0.0;
  graph.DrawSolution(0,time);
}

void ProcessamentoLocal(TPZGeoMesh &gmesh,ostream &out) {

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
	  TPZGeoEl1d  *el1d=0;
	  TPZGeoElT2d *elt2d=0;
	  TPZGeoElQ2d *elq2d=0;
	  TPZGeoElT3d *elt3d=0;
	  TPZGeoElPi3d *elpi3d=0;
	  TPZGeoElPr3d *elpr3d=0;
	  TPZGeoElC3d  *elc3d=0;
	  if(elemtype==1) el1d   = (TPZGeoEl1d    *) gel;
	  if(elemtype==2) elt2d  = (TPZGeoElT2d   *) gel;
	  if(elemtype==3) elq2d  = (TPZGeoElQ2d   *) gel;
	  if(elemtype==4) elt3d  = (TPZGeoElT3d   *) gel;
	  if(elemtype==5) elpi3d = (TPZGeoElPi3d  *) gel;
	  if(elemtype==6) elpr3d = (TPZGeoElPr3d  *) gel;
	  if(elemtype==7) elc3d  = (TPZGeoElC3d   *) gel;
	  out << "Elemento " << el->Reference()->Id() << endl;;
	  TPZManVector<REAL> sol(1);
	  TPZVec<REAL> csi(3,0.),x(3);
	  ifstream in("pontos.in");
	  in >> numpoints;
	  for(int p=0;p<numpoints;p++) {
	    in >> csi[0] >> csi[1] >>  csi[2];
	    gel->X(csi,x);
	    int var = 4;
	    if(elemtype==1) el1d->Reference()->Solution(csi,var,sol);
	    if(elemtype==2) elt2d->Reference()->Solution(csi,var,sol);
	    if(elemtype==3) elq2d->Reference()->Solution(csi,var,sol);
	    if(elemtype==4) elt3d->Reference()->Solution(csi,var,sol);
	    if(elemtype==5) elpi3d->Reference()->Solution(csi,var,sol);
	    if(elemtype==6) elpr3d->Reference()->Solution(csi,var,sol);
	    if(elemtype==7) elc3d->Reference()->Solution(csi,var,sol);
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
  cout << "\nEntre nivel inicial da malha -> 0";
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
  //retorna o n�el do elemento gel
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

REAL SolutionError(TPZCompMesh *cmesh){

  REAL Energy_error = 0.0,
    error1 = 0.0,
    error2 = 0.0;

  //erro da aproximacao FEM
  TPZManVector<REAL> flux(9);
  cmesh->EvaluateError(Solution,Energy_error,error1,error2);
  cout << "\nError Energia    : " << Energy_error << "\n\n";
  return Energy_error;
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
  LoadSolution(eixos,x,u,deriv);
  for(int i=0;i<6;i++) result[i] = u(i,0);
}

void LoadSolution(TPZFMatrix &axes,TPZVec<REAL> &X,TPZFMatrix &u,TPZFMatrix &du){

  REAL k,eps,div;
  static int key = 1;
  static REAL val;

  k = 0.0;//0.00001;
  eps = 0.4;//-1.0 �interessante
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
  u(3,0) = -( 2.*expy*(expx/div+k)*(y-1.) ) / div / eps2;//x
  u(4,0) =  ( 2.*expx*(expy/div+k)*(x-1.) ) / div / eps2;//y
  u(5,0) =  0.0;//z

  //du exact
  TPZFMatrix dur(2,6);
  dur(0,0)  =  0.0;//du/dx
  dur(1,0)  =  0.0;//du/dy

  dur(0,1)  =  0.0;//dv/dx
  dur(1,1)  =  0.0;//dv/dx

  dur(0,2)  =  -exmeydivmk*(x-1.);//dw/dx
  dur(1,2)  =  -eymexdivmk*(y-1.);//dw/dy

  dur(0,3)  =  4.*exy*(x-1.)*(y-1.) / (div*div) / (eps2*eps2);//dx/dx
  dur(1,3)  = -eymexdivmk + 2.*eymexdivmk*(y-1.)*(y-1.) / eps2;//dx/dy

  dur(0,4)  =  exmeydivmk - 2.*exmeydivmk*(x-1.)*(x-1.) / eps2;//dy/dx
  dur(1,4)  = -4.*exy*(x-1.)*(y-1.) / (div*div) / (eps2*eps2);//dy/dy

  dur(0,5)  = 0.0;//dz/dx
  dur(1,5)  = 0.0;//dz/dy

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

void LoadSolution2(TPZFMatrix &axes,TPZVec<REAL> &x,TPZFMatrix &u,TPZFMatrix &du){

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
void CoutTime(clock_t &start,char *title);
void AdaptativeProcedure(REAL error,int numiter,int marcha,int resolution,TPZMaterial *mat,TPZAnalysis &an) {
  cout << "PZAnalysis::IterativeProcessTest beginning of the iterative process, general time 0\n";
  int iter = 0,draw=0;
  TPZCompMesh *cmesh = an.Mesh();
  int numeq = cmesh->NEquations();
  cout << "\nmain::AdaptativeProcedure numero de equacoes    : " << numeq;
  int numel = NumElComp(*cmesh);
  fBegin = clock();
  an.Solution().Zero();
  an.Run();
  SolutionError(cmesh);
  cout << "TPZIterativeAnalysis::AdaptativeProcedure iteracao = " << ++iter << endl;
  CoutTime(fBegin,"TPZIterativeAnalysis::AdaptativeProcedure End system solution first iteration");
  an.LoadSolution();
  REAL normsol = Norm(an.Solution());
  REAL tolmax = 1.e15;
  REAL tol = error / sqrt(REAL(numel));
  fBegin = clock();
  int refined = AdaptiveMesh(*cmesh,tol,iter);
  CoutTime(fBegin,"TPZIterativeAnalysis::AdaptativeProcedure End adaptative procedure");

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
    CoutTime(fBegin,"TPZIterativeAnalysis::AdaptativeProcedure End system solution actual iteration");
    CoutTime(fInit,"TPZIterativeAnalysis::AdaptativeProcedure accumulated time");
    SolutionError(cmesh);
    tol = error / sqrt(REAL(numel));
    fBegin = clock();
    refined = AdaptiveMesh(*cmesh,tol,iter);
    CoutTime(fBegin,"TPZIterativeAnalysis::AdaptativeProcedure End adaptative procedure");
    cout << "TPZIterativeAnalysis::AdaptativeProcedure iteracao = " << ++iter << endl;

  }

  if( normsol > tolmax ){
    cout << "\nTPZIterativeAnalysis::AdaptativeProcedure the iterative process stopped due the great norm "
	 << "of the solution, norm solution = " << normsol << endl;
  }
  an.LoadSolution();
  CoutTime(fInit,"TPZIterativeAnalysis::AdaptativeProcedure general time of iterative process");
  if(!refined && iter <= numiter){
    cout << "\nTPZIterativeAnalysis::AdaptativeProcedure \n\n"
	 <<"\t\t\t\t\t<-> TOLERANCIA ATINGIDA <->\n\n";
  }
  if(refined == -1) cout << "\n\nPZIterativeAnalysis::AdaptativeProcedure STOP FOR MAXIM DOF\n\n";
  TPZVec<char *> scalar(2),vector(0);
  scalar[0] = "Deslocz";
  scalar[1] = "POrder";
  //scalar[1] = "Mn1";
  //scalar[2] = "Mn2";
  int dim = mat->Dimension();
  TPZDXGraphMesh graph(cmesh,dim,mat,scalar,vector);
  ofstream *dxout = new ofstream("PLACA.dx");
  cout << "\nTPZIterativeAnalysis:: out file : Placa.dx\n";
  graph.SetOutFile(*dxout);
  graph.SetResolution(resolution);
  graph.DrawMesh(dim);
  REAL time = 0.1;
  graph.DrawSolution(draw++,time);
  dxout->flush();
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
  int i,ordmax = 4;
  TPZBlock *flux;
  REAL energy,error2,error3;
  TPZVec<int> sub(0);
  TPZStack<int> nindex,pindex;
  nindex.Resize(0);
  pindex.Resize(0);
  TPZCompEl *com;
  TPZInterpolatedElement *intel;
  REAL maxerror = 0.0;
  int elmaxerr,elid,order;
  for(i=0;i<nel;i++){
    com = cmesh.ElementVec()[i];
    if(!com) continue;
    //if(com->Type() == 16) continue
    if(com->Material()->Id() < 0) continue;
    intel = dynamic_cast<TPZInterpolatedElement *>(com);
    intel->EvaluateError(Solution,energy,error2,flux,error3);
    elid = com->Reference()->Id();
    erros << "\nElemento id   = " << elid;
    erros << "\nNorma Energia = " << energy;
    if(maxerror < energy){
      maxerror = energy;
      elmaxerr = elid;
    }
    intel = dynamic_cast<TPZInterpolatedElement *>(com);
    order = intel->SideOrder(intel->NConnects()-1);
    if(energy > tol && nel < 51){//at�512 �nivel 5
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
    cout << "main::AdaptiveMesh Maximo erro da malha : " << maxerror;
    cout << "\nmain::AdaptiveMesh elemento de id = " << elid << endl;
  }
  int nsize = nindex.NElements();
  int psize = pindex.NElements(),neworder;
  if( maxerror > tol && !numdivide ) return -1;//erro uniforme n� atingido, n� h�adapta�o
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
    neworder = order < ordmax ? (order +1) : order;
    if(neworder > order) intel->PRefine(neworder);
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

//       cmesh->Print(MALHA);
//       cmesh->Reference()->Print(MALHA);
//       MALHA.flush();
//       MALHA.close();
