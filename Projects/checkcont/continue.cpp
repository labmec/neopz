///#include "pzmetis.h"
//#include "pztrnsform.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzfmatrix.h"
#include "pzelgc3d.h"
#include "pzbndcond.h"
#include "pztempmat.h"
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
#include "TPZConsLawTest.h"
//#include "TPZRefPattern.h"
#include "TPZCompElDisc.h"
#include "TPZShapeDisc.h"
#include "TPZInterfaceEl.h"
//#include "TPZJacobMat.h"
//#include "TPZJacobStrMatrix.h"

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <iostream>

//#define NOTDEBUG

void CriacaoDeNos(int nnodes,double lista[20][3]);
void Linha();
void Hexaedro();
void Tetraedro();
void Piramide();
void Prisma();
void ContagemDeElementos();
void CycleRefinements(TPZCompMesh &cm, int numcycles, int minel, int maxel, ofstream &out);
int IsGroup(int ind, TPZCompEl *cel, TPZVec<int> &indgroup);

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


void BCMenos1(TPZVec<REAL> &x,TPZVec<REAL> &result);
void BCMais1(TPZVec<REAL> &x,TPZVec<REAL> &result);
void CC_Poisson(TPZMaterial *mat);
void CC_Projecao(TPZMaterial *mat);
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
static int problemtype;

/* class TPZGeoElC3d; */
/* class TPZGeoElT3d; */
/* TPZGeoElC3d::SetCreateFunction(TPZCompElDisc::CreateC3Disc); */
/* TPZGeoElT3d::SetCreateFunction(TPZCompElDisc::CreateT3Disc); */
/* TPZGeoElPi3d::SetCreateFunction(TPZCompElDisc::CreatePi3Disc); */
/* TPZGeoElPr3d::SetCreateFunction(TPZCompElDisc::CreatePr3Disc); */
/* TPZGeoElT2d::SetCreateFunction(TPZCompElDisc::CreateT2Disc); */
/* TPZGeoElQ2d::SetCreateFunction(TPZCompElDisc::CreateQ2Disc); */


int main() {


  int interfdim;
  cout << "Entre com a dimensao da interface : ";
  cin >> interfdim;
  TPZCompElDisc::gInterfaceDimension = interfdim;
  int tipo;
  cout << "\nElemento [0:linha][1:tetraedro][2:piramide][3:prisma][4:hexaedro]\n\t";
  cin >> tipo;
  if(tipo==0){
    CriacaoDeNos(2,linha);// nodes
    Linha();// elementos
  }  
  if(tipo==1){
    CriacaoDeNos(4,tetraedro);// nodes
    Tetraedro();// elementos
  }
  if(tipo==2){
    CriacaoDeNos(5,piramide);// nodes
    Piramide();// elementos
  }
  if(tipo==3){
    CriacaoDeNos(6,prisma);// nodes
    Prisma();// elementos
  }
  if(tipo==4){
    CriacaoDeNos(8,hexaedro);// nodes
    Hexaedro();// elementos
  }
  
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
    ContagemDeElementos();
  }

  if(1){
    cout << "\nmain::Imprime malhas\n";
    gmesh->Print(outgm);
    cmesh->Print(outgm);
    outgm.flush();
    outgm.close();
  }  

  if(0){
    cout << "\nmain::Ajuste no contorno e imprime malhas\n";
    cmesh->AdjustBoundaryElements();
    gmesh->Print(outgm);
    cmesh->Print(outgm);
    outgm.flush();
    outgm.close();
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
      TPZStepSolver solver;
      solver.SetDirect(ELDLt);
      an.SetSolver(solver);
      an.Run(outgm);
      if(0)
	{//pos-processamento
	  TPZVec<char *> scalar(1),vector(0);
	  scalar[0] = "Solution";
	  //vector[0] = "Solution";//"Tension9";
	  char plot[] = "conslaw.dx";
	  an.DefineGraphMesh(3,scalar,vector,plot);
	  an.PostProcess(0);
	}
    }
    
    //ofstream outan("analysis.out");
    an.Print("FEM SOLUTION ",outgm);

    REAL estimate=0.,H1_error=0.,L2_error=0.;

    if(1){//erro da aproximacao FEM  
      TPZManVector<REAL> flux(9);
      cmesh->EvaluateError(Solution,H1_error,L2_error,estimate);
      outgm << "\n\n\n";    
      outgm << "L2 Error  : " << L2_error << endl;
      outgm << "Semi Norm : " << estimate << endl;
      outgm << "Energy    : " << H1_error << endl;
      outgm.flush();
    }
    cout << "L2 Error  : " << L2_error << endl;
    cout << "Semi Norm : " << estimate << endl;
    cout << "Energy    : " << H1_error << endl;

    if(0){//pos-processamento local
      PostProcess(*gmesh,outgm);
      outgm.flush();  
    }

    ContagemDeElementos();
    if(problemtype == 1) cout << "\n\n------- Problema Projecao -------\n\n";
    if(problemtype == 2) cout << "\n------- Problema Poisson -------\n\n";
  }//if(0)

  if(cmesh) delete cmesh;
  if(gmesh) delete gmesh;
  //AvisoAudioVisual();
  return 0;
}

void ContagemDeElementos(){

  int poin=0,line=0,tria=0,quad=0,tetr=0,pira=0,pris=0,hexa=0,disc=0,inte=0;
  int nelem = cmesh->ElementVec().NElements();
  int k,totel=0,bcel=0,nivel = 0,nivmax=0;
  for(k=0;k<nelem;k++){
    TPZCompEl *comp = cmesh->ElementVec()[k];
    if(!comp) continue;
    totel++;
    if(comp->Reference()->MaterialId() < 0) bcel++;
    nivel = comp->Reference()->Level();
    if(nivmax < nivel) nivmax = nivel;
/*     if(comp->Type() == 00) poin++; */
/*     if(comp->Type() == 01) line++; */
/*     if(comp->Type() == 02) tria++; */
/*     if(comp->Type() == 03) quad++; */
/*     if(comp->Type() == 04) tetr++; */
/*     if(comp->Type() == 05) pira++; */
/*     if(comp->Type() == 06) pris++; */
/*     if(comp->Type() == 07) hexa++; */
    if(comp->Type() == 16) disc++;
    if(comp->Type() == 17) inte++;
  }
  nelem = gmesh->ElementVec().NElements();
  int total=0,nivmax2=0;
  for(k=0;k<nelem;k++){
    TPZGeoEl *geo = gmesh->ElementVec()[k];
    if(!geo) continue;
    total++;
    nivel = geo->Level();
    if(nivmax2 < nivel) nivmax2 = nivel;
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
  cout << "\nNivel maximo comput. atingido      : " << nivmax;              ;
  cout << "\nNivel maximo geomet. atingido      : " << nivmax << endl << endl;
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

void Linha(){
 
  CriacaoDeNos(2,linha);
  //elemento de volume
  TPZVec<int> nodes;
  nodes.Resize(2);
  nodes[0] = 0;
  nodes[1] = 1;
  TPZGeoEl1d *elg1d = new TPZGeoEl1d(nodes,1,*gmesh);
  TPZGeoEl1d::SetCreateFunction(TPZCompElDisc::Create1dDisc);//de volume
  TPZGeoElPoint::SetCreateFunction(TPZCompElDisc::CreatePointDisc);//não de interface
  gmesh->BuildConnectivity();
  TPZVec<REAL> B(1);
  B[0] = 1.0;
  int nummat = 1;
  int dif_artif = 1;
  REAL delta_t = 0.1;
  int dim = 1;
  TPZMaterial *mat = new TPZConsLawTest(nummat,B,dif_artif,delta_t,dim);
  //mat->SetForcingFunction(Function);
  cmesh->InsertMaterialObject(mat);
  //condições de contorno  
  TPZBndCond *bc;
  REAL big = 1.e12;
  TPZFMatrix val1(1,1,0.),val2(1,1,0.);
  
  //CC no canto -1
  val1.Zero();//para uso CC tipo 2
  val2.Zero();
  TPZGeoElBC(elg1d,0,-1,*gmesh);
  bc = mat->CreateBC(-1,0,val1,val2);
  bc->SetForcingFunction(BCMenos1);
  cmesh->InsertMaterialObject(bc);
  
  //CC no canto 1
  TPZGeoElBC(elg1d,1,-2,*gmesh);
  bc = mat->CreateBC(-2,0,val1,val2);
  bc->SetForcingFunction(BCMais1);
  cmesh->InsertMaterialObject(bc); 

  int ordem=1;
  cout << "\nEntre ordem -> 1";
  //cin >> ordem;
  TPZCompEl::gOrder = ordem;
  cout << endl;
  cmesh->AutoBuild();
}

void BCMenos1(TPZVec<REAL> &x,TPZVec<REAL> &result){

    result[0] = 1.0;
}

void BCMais1(TPZVec<REAL> &x,TPZVec<REAL> &result){

    result[0] =  2.0;
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
  TPZGeoElPi3d *elg3d = new TPZGeoElPi3d(nodes,100,*gmesh);
  //elementos de contorno (não de interface)
  nodes.Resize(4);//o primeiro é quadrilátero
  int NODES[5][4] = {{0,1,2,3},{0,1,4,-1},{1,2,4,-1},{3,2,4,-1},{0,3,4,-1}};
  int i;
  for(i=0;i<5;i++){
    nodes[0] = NODES[i][0];
    nodes[1] = NODES[i][1];
    nodes[2] = NODES[i][2];
    if(i==0){
      nodes[3] = NODES[i][3];
      TPZGeoElQ2d *elg2d = new TPZGeoElQ2d(nodes,100,*gmesh);
      nodes.Resize(3);//os próximos são triângulos
      continue;
    }
    TPZGeoElT2d *elg2d = new TPZGeoElT2d(nodes,100,*gmesh);
  }
  TPZGeoElPi3d::SetCreateFunction(TPZCompElDisc::CreatePi3Disc);
  TPZGeoElT3d::SetCreateFunction(TPZCompElDisc::CreateT3Disc);
  TPZGeoElT2d::SetCreateFunction(TPZCompElDisc::CreateT2Disc);
  TPZGeoElQ2d::SetCreateFunction(TPZCompElDisc::CreateQ2Disc);
  gmesh->BuildConnectivity();
  TPZMaterial *mat = new TPZMatPoisson3d(100,3);
  mat->SetForcingFunction(Function);
  cmesh->InsertMaterialObject(mat);
  //condições de contorno  
  //if(problemtype == 1) CC_Projecao(mat);
  //if(problemtype == 2) CC_Poisson(mat);
  int ordem=2;
  //cout << "\nEntre ordem -> ";
  //cin >> ordem;
  TPZCompEl::gOrder = ordem;
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
  TPZGeoElPr3d *elg3d = new TPZGeoElPr3d(nodes,100,*gmesh);
  //elementos de contorno (não de interface)
  nodes.Resize(3);//o primeiro é triângulo
  int NODES[5][4] = {{0,1,2,-1},{0,1,4,3},{1,2,5,4},{0,2,5,3},{3,4,5,-1}};
  int i;
  for(i=0;i<5;i++){
    nodes.Resize(3);
    nodes[0] = NODES[i][0];
    nodes[1] = NODES[i][1];
    nodes[2] = NODES[i][2];
    if(i>0 && i<4){
      nodes.Resize(4);
      nodes[3] = NODES[i][3];
      TPZGeoElQ2d *elg2d = new TPZGeoElQ2d(nodes,100,*gmesh);
      continue;
    }
    TPZGeoElT2d *elg2d = new TPZGeoElT2d(nodes,100,*gmesh);
  }
  TPZGeoElPr3d::SetCreateFunction(TPZCompElDisc::CreatePr3Disc);
  TPZGeoElT2d::SetCreateFunction(TPZCompElDisc::CreateT2Disc);
  TPZGeoElQ2d::SetCreateFunction(TPZCompElDisc::CreateQ2Disc);
  gmesh->BuildConnectivity();
  TPZMaterial *mat = new TPZMatPoisson3d(100,3);
  mat->SetForcingFunction(Function);
  cmesh->InsertMaterialObject(mat);
  //condições de contorno  
  //if(problemtype == 1) CC_Projecao(mat);
  //if(problemtype == 2) CC_Poisson(mat);
  int ordem=2;
  //cout << "\nEntre ordem -> ";
  //cin >> ordem;
  TPZCompEl::gOrder = ordem;
  cout << endl;
  cmesh->AutoBuild();
}

void Tetraedro(){
 
  //elemento de volume
  TPZVec<int> nodes;
  nodes.Resize(4);
  nodes[0] = 0;
  nodes[1] = 1;
  nodes[2] = 2;
  nodes[3] = 3;
  TPZGeoElT3d *elg3d = new TPZGeoElT3d(nodes,100,*gmesh);
  //elementos de contorno (não de interface nem BC)
  nodes.Resize(3);
  int NODES[4][3] = {{0,1,2},{0,1,3},{1,2,3},{0,2,3}};
  int i;
  for(i=0;i<4;i++){
    nodes[0] = NODES[i][0];
    nodes[1] = NODES[i][1];
    nodes[2] = NODES[i][2];
    TPZGeoElT2d *elg2d = new TPZGeoElT2d(nodes,100,*gmesh);
  }
  //construtor descontínuo
  TPZGeoElT3d::SetCreateFunction(TPZCompElDisc::CreateT3Disc);
  TPZGeoElPi3d::SetCreateFunction(TPZCompElDisc::CreatePi3Disc);
  TPZGeoElQ2d::SetCreateFunction(TPZCompElDisc::CreateQ2Disc);
  TPZGeoElT2d::SetCreateFunction(TPZCompElDisc::CreateT2Disc);
  gmesh->BuildConnectivity();
  TPZMaterial *mat = new TPZMatPoisson3d(100,3);
  mat->SetForcingFunction(Function);
  cmesh->InsertMaterialObject(mat);
  int ordem=2;
  //cout << "\nEntre ordem -> ";
  //cin >> ordem;
  TPZCompEl::gOrder = ordem;
  cout << endl;
  cmesh->AutoBuild();
}

void Hexaedro(){
 
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
  TPZGeoElC3d *elg3d = new TPZGeoElC3d(nodes,100,*gmesh);
  //elementos de contorno (não de interface nem BC)
  nodes.Resize(4);
  int NODES[6][4] = {{0,1,2,3},{0,1,5,4},{1,2,6,5},{3,2,6,7},{0,3,7,4},{4,5,6,7}};
  int i;
  for(i=0;i<6;i++){
    nodes[0] = NODES[i][0];
    nodes[1] = NODES[i][1];
    nodes[2] = NODES[i][2];
    nodes[3] = NODES[i][3];
    TPZGeoElQ2d *elg2d = new TPZGeoElQ2d(nodes,100,*gmesh);
  }
  //construtor atual
  //TPZCompEl *(*fpq)(TPZGeoElQ2d *geoel,TPZCompMesh &mesh,int &index);
  //TPZCompEl *(*fpt)(TPZGeoElT2d *geoel,TPZCompMesh &mesh,int &index);
  //fpq = TPZGeoElQ2d::fp;
  //fpt = TPZGeoElT2d::fp;
  //construtor descontínuo
  TPZGeoElC3d::SetCreateFunction(TPZCompElDisc::CreateC3Disc);
  TPZGeoElQ2d::SetCreateFunction(TPZCompElDisc::CreateQ2Disc);
  gmesh->BuildConnectivity();
  TPZMaterial *mat = new TPZMatPoisson3d(100,3);
  mat->SetForcingFunction(Function);
  cmesh->InsertMaterialObject(mat);
  int ordem=2;
  cout << "\nEntre ordem -> ";
  //cin >> ordem;
  TPZCompEl::gOrder = ordem;
  cout << endl;
  cmesh->AutoBuild();
  //TPZGeoElQ2d::SetCreateFunction(fpq);//construtores
  //TPZGeoElT2d::SetCreateFunction(fpt);//contínuos usuais
}

void CC_Projecao(TPZMaterial *mat){

  cout << "main nao usar\n";

   TPZBndCond *bc;
   REAL big = 1.e12;
   TPZFMatrix val1(1,1,0.),val2(1,1,0.);       
   TPZGeoEl *elg0 = gmesh->ElementVec()[0];

   //CC na face 20
   val1.Zero();//para uso CC tipo 2
   val2.Zero();
   TPZGeoElBC(elg0,20,-1,*gmesh);
   bc = mat->CreateBC(-1,0,val1,val2);
   bc->SetForcingFunction(Function);
   cmesh->InsertMaterialObject(bc);

   //CC na face 21
   TPZGeoElBC(elg0,21,-2,*gmesh);
   bc = mat->CreateBC(-2,0,val1,val2);
   bc->SetForcingFunction(Function);
   cmesh->InsertMaterialObject(bc);

   //CC na face 22
   TPZGeoElBC(elg0,22,-3,*gmesh);
   bc = mat->CreateBC(-3,0,val1,val2);
   bc->SetForcingFunction(Function);
   cmesh->InsertMaterialObject(bc);

   //CC na face 23
   TPZGeoElBC(elg0,23,-4,*gmesh);
   bc = mat->CreateBC(-4,0,val1,val2);
   bc->SetForcingFunction(Function);
   cmesh->InsertMaterialObject(bc);

   //CC na face 24
   TPZGeoElBC(elg0,24,-5,*gmesh);
   bc = mat->CreateBC(-5,0,val1,val2);
   bc->SetForcingFunction(Function);
   cmesh->InsertMaterialObject(bc);

   //CC na face 25
   TPZGeoElBC(elg0,25,-6,*gmesh);
   bc = mat->CreateBC(-6,0,val1,val2);
   bc->SetForcingFunction(Function);
   cmesh->InsertMaterialObject(bc);
}

void CC_Poisson(TPZMaterial *mat){

  cout << "main nao usar\n";

  //nulo nos lados contidos nos eixos
  TPZBndCond *bc;
  REAL big = 1.e12;
  TPZFMatrix val1(1,1,0.),val2(1,1,0.);       
  TPZGeoEl *elg0 = gmesh->ElementVec()[0];
  
  //CC na face 20
  val1.Zero();//para uso CC tipo 2
  val2.Zero();
  TPZGeoElBC(elg0,20,-1,*gmesh);
  bc = mat->CreateBC(-1,0,val1,val2);
  bc->SetForcingFunction(BCPoisson);
  cmesh->InsertMaterialObject(bc);
  
  //CC na face 21
  TPZGeoElBC(elg0,21,-2,*gmesh);
  bc = mat->CreateBC(-2,0,val1,val2);
  bc->SetForcingFunction(BCPoisson);
  cmesh->InsertMaterialObject(bc);
  
  //CC na face 22
  TPZGeoElBC(elg0,22,-3,*gmesh);
  bc = mat->CreateBC(-3,1,val1,val2);
  bc->SetForcingFunction(BCPoisson);
  cmesh->InsertMaterialObject(bc);
  
  //CC na face 23
  TPZGeoElBC(elg0,23,-4,*gmesh);
  bc = mat->CreateBC(-4,1,val1,val2);
  bc->SetForcingFunction(BCPoisson);
  cmesh->InsertMaterialObject(bc);

  //CC na face 24
  TPZGeoElBC(elg0,24,-5,*gmesh);
  bc = mat->CreateBC(-5,0,val1,val2);
  bc->SetForcingFunction(BCPoisson);
  cmesh->InsertMaterialObject(bc);

  //CC na face 25
  TPZGeoElBC(elg0,25,-6,*gmesh);
  bc = mat->CreateBC(-6,1,val1,val2);
  bc->SetForcingFunction(BCPoisson);
  cmesh->InsertMaterialObject(bc);
}

void BCPoisson(TPZVec<REAL> &x,TPZVec<REAL> &result){

  //solução u = x2+y2+z2, domínio cubo unitário
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
  //no plano Xi=0 a derivada normal depende de Xi, logo é nula
  //a condição imposta é de Dirichlet, a CC de Neumman é nula
  if( fabs(x[0]) < 1.e-8 || fabs(x[1]) < 1.e-8 || fabs(x[2]) < 1.e-8)
    result[0] = x[0]*x[0]+x[1]*x[1]+x[2]*x[2];
}

void Function(TPZVec<REAL> &x,TPZVec<REAL> &result){
  //sol u = x*x + y*y + z*z , -laplaciano(u) = f  
  if(problemtype == 1) result[0]  = x[0]*x[0]+x[1]*x[1]+x[2]*x[2];
  if(problemtype == 2) result[0] = -6.0;//Laplaciano
}

void Solution(TPZVec<REAL> &x,TPZVec<REAL> &result,TPZFMatrix &deriv){

  result[0]  = 0.;
  deriv(0,0) = 0.;
  deriv(1,0) = 0.;
  deriv(2,0) = 0.;
}

void Divisao (TPZCompMesh *cmesh,int key){
  
  int k=0;
  if(key < 0) return;
  TPZVec<int> csub(0);
  int n1=1,n2;
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
    if(cpel && el < nelc && cpel->Type() == 17){
      PZError << "main::Divisao elemento interface (nao foi dividido!)\n\n";
      cout << "Elementos divissiveis:\n";
      for(el=0;el<nelc;el++) {
	cpel = cmesh->ElementVec()[el];
	if(cpel && cpel->Type() != 17){
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


static REAL point[4][3] = {{0.,0.,0.},{.1,.2,.3},{.1,.1,.1},{-.1,-.1,-.1}};
void PostProcess(TPZGeoMesh &gmesh,ostream &out) {

  int numpoints;
  int nel = gmesh.Reference()->ElementVec().NElements();
  int count = -1;
  while(count++<nel){
    for(int iel=0;iel<nel;iel++) {
      if(!gmesh.Reference()->ElementVec()[iel]) continue;      
      int elemtype = gmesh.Reference()->ElementVec()[iel]->Type();
      if(elemtype==0) continue;
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

	  for(int p=0;p<4;p++) {
	    csi[0] = point[p][0]; 
	    csi[1] = point[p][1]; 
	    csi[2] = point[p][2]; 
	    gel->X(csi,x);
	    if(elemtype==1) el1d->Reference()->Solution(csi,0,sol);
	    if(elemtype==2) elt2d->Reference()->Solution(csi,0,sol);
	    if(elemtype==3) elq2d->Reference()->Solution(csi,0,sol);
	    if(elemtype==4) elt3d->Reference()->Solution(csi,0,sol);
	    if(elemtype==5) elpi3d->Reference()->Solution(csi,0,sol);
	    if(elemtype==6) elpr3d->Reference()->Solution(csi,0,sol);
	    if(elemtype==7) elc3d->Reference()->Solution(csi,0,sol);
	    out << "solucao em x    = " << x[0] << ' ' << x[1] << ' ' << x[2] << endl;
	    out << "               u = " << sol[0] << endl;	    
	  }
	} else {
	  continue;
	}
      }
    }
  }
}

void TestShapeDescontinous(){

  int nel = cmesh->ElementVec().NElements(),i;
  REAL C[3];
  if(0){
    for(i=0;i<nel;i++){
      TPZCompEl *comp = (cmesh->ElementVec()[i]);
      TPZCompElDisc *disc = dynamic_cast<TPZCompElDisc *> (comp);
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
  int const N=1;
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

  int maxBC=0;
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
  int elbc=0;
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
	if(cel->Type() == 17) continue;
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
	  realdiv++;//todos os divididos que não eram agrupados
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
	  cm.CoarsenDisc(subindex,elindex);
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
