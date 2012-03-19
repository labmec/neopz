/**
 * @file
 * @brief Contains a tutorial example using a discontinuous approximation function.
 * @note c = sqrt(gama*p/ro)  velocidade do som
*/

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
#include "TPZGeoElement.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzfmatrix.h"
#include "pzbndcond.h"

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
#include "pzelmat.h"
#include "pzelasmat.h"
#include "pzmattest.h"
#include "pzmat1dlin.h"
#include "pzmat2dlin.h"
#include "pzpoisson3d.h"
#include "pzmaterial.h"
#include "pzconslaw.h"
#include "pzeulerconslaw.h"
#include "TPZDiffusionConsLaw.h"
#include "TPZConsLawTest.h"
#include "TPZCompElDisc.h"
#include "TPZShapeDisc.h"
#include "TPZInterfaceEl.h"
#include "TPZExtendGridDimension.h"
#include "pzreal.h"
#include "TPZNLMultGridAnalysis.h"
#include "pzflowcmesh.h"
#include "pzdxmesh.h"
#include "pzmganalysis.h"

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

#include <iostream>
#include <ostream>
#include <cmath>

using namespace std;

static double quadrado[4][3] = { {0.,0.,0.},{4.,0.,0.},{4.,4.,0.},{0.,4.,0.} };

static REAL p1=0.7,p2=0.902026,p3=1.80405,p4=2.96598,p5=3.2,p6=4.12791,p7=0.22;

static double novequads[15][3] = { {0.,0.,0.},{p1+p1/4.,0.,0.},{p3,0.,0.},
				   {p5-p7/1.5,0.,0.},{p6,0.,0.},{0.,p7,0.},{p1,p7,0.},
				   {p2,.5,0.},{p3,.75,0.},{p4,.5,0.},{p5,p7,0.},
				   {p6,p7,0.},{0.,1.,0.},{p3,1.,0.},{p6,1.,0.}};

static double novecubos[30][3] = { {0.,0.,0.},{p1+p1/4.,0.,0.},{p3,0.,0.},
				   {p5-p7/1.5,0.,0.},{p6,0.,0.},{0.,p7,0.},
				   {p1,p7,0.},{p2,.5,0.},{p3,.75,0.},{p4,.5,0.},
				   {p5,p7,0.},{p6,p7,0.},{0.,1.,0.},{p3,1.,0.},
				   {p6,1.,0.},{0.,0.,0.3},{p1+p1/4.,0.,0.3},
				   {p3,0.,0.3},{p5-p7/1.5,0.,0.3},{p6,0.,0.3},
				   {0.,p7,0.3},{p1,p7,0.3},{p2,.5,0.3},{p3,.75,0.3},
				   {p4,.5,0.3},{p5,p7,0.3},{p6,p7,0.3},
				   {0.,1.,0.3},{p3,1.,0.3},{p6,1.,0.3} };

static double quadrilatero2[5][3] = { {0.,0.,0.},{1.80405,0.,0.},{4.12791,0.,0.},{0.,1.,0.},{4.12791,1.,0.} };

void AgrupaList(TPZVec<int> &accumlist,int nivel,int &numaggl);
void SetDeltaTime(TPZMaterial *mat,REAL deltaT);
void CriacaoDeNos(int nnodes,double lista[20][3]);
void SetDeltaTime(TPZMaterial *mat,TPZCompMesh *cmesh);
TPZMaterial *NoveCubos(int grau);
TPZMaterial *NoveQuadrilateros(int grau);
TPZMaterial *ProblemaQ2D1El(int grau);
TPZMaterial *ProblemaT2D(int grau);
TPZMaterial *TresTriangulos(int grau);
TPZMaterial *Triangulo(int grau);
TPZMaterial *Quadrilatero(int grau);
TPZMaterial *Hexaedro(int grau);
TPZMaterial *TresPrismas(int grau);
TPZMaterial *FluxConst3D(int grau);
TPZMaterial *FluxConst2D(int grau);
void ContagemDeElementos(TPZMaterial *mat);
void FileNB(TPZGeoMesh &gmesh,std::ostream &out,int var);
void Function(TPZVec<REAL> &x,TPZVec<REAL> &result);
void PostProcess(TPZCompMesh &cmesh,std::ostream &out,int var);
void Divisao (TPZCompMesh *cmesh);
void NivelDivide(TPZCompMesh *cmesh);
void SequenceDivide2();
void SequenceDivide(int fat[100],int numbel);
void TestShapesDescontinous();
static clock_t end;
void CoutTime(clock_t &start);

static TPZGeoMesh *gmesh = new TPZGeoMesh;
static TPZCompMesh *cmesh = new TPZFlowCompMesh(gmesh),
*aggcmesh = new TPZFlowCompMesh(gmesh);
static int grau = 0;
static int tipo;
static int problem=0;;
static REAL CFL=-1.0,deltaT = -1.0;
static REAL gama = 1.4;
static int meshdim = -1;


int main() {

//   TPZMGAnalysis::main();
//   return 0;

  gmesh->SetName("\n\t\t\t* * * INITIAL GEOMETRIC MESH * * *\n\n");
  cmesh->SetName("\n\t\t\t* * * INITIAL COMPUTATIONAL MESH * * *\n\n");
  ofstream outgm("mesh.out");

  cout << "\nChoice\n"
       << "\t[0: ThreeTriangules]\n"
       << "\t[5: NineQuadrilateral]\n"
       << "\t[6: NineCubes]\n"
       << "\t\t\t";

  //cin >> tipo;
  tipo = 0;
  problem = tipo;
  cout << "\nDegree of the interpolation space, p = 0,1,2,3,... ";
  cin >> grau;
  //grau = 0;
  TPZCompEl::SetgOrder(grau);
  TPZConservationLaw *mat;
  if(tipo==0){
    mat = dynamic_cast<TPZConservationLaw *>(TresTriangulos(grau));
    meshdim = 2;
    problem = 0;
  } else
  if(tipo==5){
    mat = dynamic_cast<TPZConservationLaw *>(NoveQuadrilateros(grau));
    meshdim = 2;
    problem = 0;
  } else if(tipo==6){
    mat = dynamic_cast<TPZConservationLaw *>(NoveCubos(grau));
    meshdim = 3;
  } else {
    //INCOMPLETO s� interesa o dom�nio quadrado
    mat = dynamic_cast<TPZConservationLaw *>(ProblemaQ2D1El(grau));
    meshdim = 2;
    problem = -1;
  }
  if(1){
    cout << "\neuler.c::main verificando a consistencia da malha de interfaces\t";
    if(TPZInterfaceElement::main(*cmesh)){
      cout << "->\tOK!";
    } else {
      cout << "->\tPROBLEMAS COM INTERFACES\n\n";
      ContagemDeElementos(mat);
      return 0;
    }
  }
  if(0){
    cout << "\nmain::Imprime malhas\n";
    gmesh->Print(outgm);
    cmesh->Print(outgm);
    outgm.flush();
  }

  if(0) NivelDivide(cmesh);

  cout << "\nmain:: entre CFL (si nulo sera calculado) -> \n";
  cin >> CFL;//CFL = 0.3;
  TPZDiffusionConsLaw::fCFL = CFL;
//   SetDeltaTime(mat,cmesh);
//   REAL timestep = mat->TimeStep();
//   cout << "\nmain::Passo inicial de tempo -> " << timestep << endl;
//   mat->SetTimeStep(-1.0);//para reiniciar com SetDeltaTime(mat,cmesh);
  cout << "main:: delta (si nulo sera calculado) -> \n";
  cin >> TPZDiffusionConsLaw::fDelta;
  if(0){
    cout << "\nmain::Ajuste no contorno e imprime malhas\n";
    cmesh->AdjustBoundaryElements();
    if(1){
      gmesh->Print(outgm);
      cmesh->Print(outgm);
      outgm.flush();
    }
  }

  //com matriz n�o sim�trica e ELU 2D e 3D convergen
  int numat = mat->Id(),alg;
  //ofstream outan("PostProcess.out");
  cout << "\nmain:: algoritmo a uma ou duas malhas?\n";
  cin >> alg;
  if(1){
    cmesh->SetDimModel(meshdim);
    TPZNonLinMultGridAnalysis mgnlan(cmesh);
    mgnlan.SetAnalysisFunction(SetDeltaTime);
    if(alg == 1){
      mgnlan.OneGridAlgorithm(outgm,numat);
    } else if(alg == 2){
      mgnlan.TwoGridAlgorithm(outgm,numat);
    }
    int nm = mgnlan.NMeshes();
    TPZCompMesh *mesh = mgnlan.IMesh(nm-1);//malha fina
    if(1){
      mesh->Print(outgm);
      mesh->Reference()->Print(outgm);
      outgm.flush();
    }
    ContagemDeElementos(mat);
    //PostProcess(*mgnlan.IMesh(1),outan,0);//malha grossa
    //PostProcess(*mgnlan.IMesh(2),outan,0);//malha fina
  }
  //outan.flush();
  //outan.close();
  outgm.close();
  //if(cmesh) delete cmesh;//foi apagada no multigrid
  if(gmesh) delete gmesh;
  return 0;
}

//* FIM_MAIN * FIM_MAIN * FIM_MAIN * FIM_MAIN * FIM_MAIN * FIM_MAIN * FIM_MAIN
//* FIM_MAIN * FIM_MAIN * FIM_MAIN * FIM_MAIN * FIM_MAIN * FIM_MAIN * FIM_MAIN

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

void SetDeltaTime(TPZMaterial *mat,REAL deltaT){


  TPZConservationLaw *law = dynamic_cast<TPZConservationLaw *>(mat);

  law->SetTimeStep(deltaT);

}

void Divisao (TPZCompMesh *cmesh){

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
    if(cpel && el < nelc && cpel->Type() == EInterface){
      PZError << "main::Divisao elemento interface (nao foi dividido!)\n\n";
      cout << "Elementos divissiveis:\n";
      for(el=0;el<nelc;el++) {
	cpel = cmesh->ElementVec()[el];
	if(cpel && cpel->Type() != EInterface){
	  TPZGeoEl *gel = cpel->Reference();
	  if(gel) cout << gel->Id() << ",";
	}
      }
    } else {
      if(!el || el < nelc) cmesh->Divide(el,csub,0);
      else {
	cout << "main::Divisao elemento sem referencia\n";
	ContagemDeElementos(0);
      }
      n1 = 1;
    }
  }
}

using namespace pzshape;

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
  TPZFMatrix<REAL> phi,dphi;
  //  int const N=1;
  TPZShapeDisc::Shape1D(C[0],X0,X,degree,phi,dphi);
  phi.Print("Uni-dimensional",cout);
  dphi.Print("Uni-dimensional",cout);
  phi.Resize(0,0);
  dphi.Resize(0,0);
  TPZShapeDisc::Shape2D(C[1],X0,X,degree,phi,dphi,TPZShapeDisc::ETensorial);
  phi.Print("Bi-dimensional",cout);
  dphi.Print("Bi-dimensional",cout);
  phi.Resize(0,0);
  dphi.Resize(0,0);
  TPZShapeDisc::Shape3D(C[2],X0,X,degree,phi,dphi,TPZShapeDisc::ETensorial);
  phi.Print("Tri-dimensional",cout);
  dphi.Print("Tri-dimensional",cout);
}

void ContagemDeElementos(TPZMaterial *mat){

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
    if(comp->Type() == EDiscontinuous) disc++;
    if(comp->Type() == EInterface) inte++;
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
  cout << "\nGrau do espa�o de interpola��o     : " << grau;
  cout << "\nNivel maximo comput. atingido      : " << nivmax;
  cout << "\nNivel maximo geomet. atingido      : " << nivmax << endl << endl;
  if(mat){
    cout << "\nPropriedades materiais             : ";
    mat->Print();
    cout << "Time Step = : " << deltaT << endl;
  }
  cout << "\nDeltaX                             : " <<  cmesh->DeltaX() << endl;
  cout << "\nLesserEdgeOfEl                     : " <<  cmesh->LesserEdgeOfMesh() << endl;
  cout << "\nMaximumRadiusOfEl                  : " <<  cmesh->MaximumRadiusOfMesh() << endl;
}

int Nivel(TPZGeoEl *gel);
void NivelDivide(TPZCompMesh *cmesh){

  TPZVec<int> csub(0);
  int nivel;
  cout << "\nmain::Divisao todos os elementos da malha serao divididos!\n";
  cout << "\nmain::Divisao Nivel da malha final ? : ";
  cin >> nivel;
  cout << "\nNivel da malha a ser atingido = " << nivel << endl;
  int nelc = cmesh->ElementVec().NElements();
  int el,actual;
  TPZCompEl *cpel;
  TPZGeoEl *gel;
  el = -1;
  while(++el<nelc) {
    cpel = cmesh->ElementVec()[el];
    if(!cpel) continue;
    if(cpel->Type() == EInterface) continue;
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
			   {-0.8,-0.8,00.8},{0.8,-0.8,00.8},{0.8,0.8,00.8},{-0.8,0.8,00.8},
			   {0.,0.,0.} };//hexaedro
void PostProcess(TPZCompMesh &cmesh,std::ostream &out,int var) {

  int nel = cmesh.ElementVec().NElements();
  TPZGeoMesh *gmesh = cmesh.Reference();
  if(nel > 1000){
    cout << "main::PostProcess mas de 10000 elementos -> processa 2000\n";
  }
  int idmax = 0,dim=cmesh.Dimension(),finish=-1,i;
  for(int iel=0;iel<nel;iel++){//procurando o id mais alto da lista
    if(++finish >= 2000) return;
    TPZCompEl *cel = cmesh.ElementVec()[iel];
    if(!cel) continue;
    TPZGeoEl *el = gmesh->ElementVec()[iel];
    int id = el->Id();
    if(id > idmax) idmax = id;
  }
  for(int iel=0;iel<nel;iel++) {
    TPZCompEl *el = cmesh.ElementVec()[iel];
    if(!el) continue;
    if(el->Dimension() != dim) continue;
//    int elemtype = el->Type();
    TPZGeoEl *gel = el->Reference();
    if(el && gel) {
//      int nsides = gel->NSides();
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
        gel->Reference()->Solution(csi,var,sol);
	int size = sol.NElements();
	out << "solucao em x    = " << x[0] << ' ' << x[1] << ' ' << x[2];
	out << " -> u = ";
	for(i=0;i<size;i++)
	  out << sol[i] << " ";
	out << endl;
      }
    }
  }
}

void Ordena(TPZVec<REAL> &coordx,TPZVec<int> &sort);
void FileNB(TPZGeoMesh &gmesh,std::ostream &out,int var) {

  //  int numpoints;
  int nel = gmesh.Reference()->ElementVec().NElements();
  int idmax = 0,dim,chega=1,finish=-1;
  for(int iel=0;iel<nel;iel++){//procurando o id mais alto da lista
    if(++finish >= 5000) return;
    TPZCompEl *cel = gmesh.Reference()->ElementVec()[iel];
    if(!cel) continue;
    TPZGeoEl *el = gmesh.ElementVec()[iel];
    if(chega && cel->Type()==EDiscontinuous) {dim = el->Dimension(); chega = 0;}
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
      if(elemtype==EInterface) continue;//interface
      TPZCompEl *el = gmesh.Reference()->ElementVec()[iel];
      if(el->Material()->Id() < 0) continue;
      //s� elementos de volume
      TPZGeoEl *gel = el->Reference();
      if(el && gel) {
	if(gel->Id()==count){
	  //int nsides = gel->NSides();
	  TPZManVector<REAL> sol(1);
	  TPZVec<REAL> csi(3,0.),x(3);
	  //int var = 1;//densidade
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
            gel->Reference()->Solution(csi,var,sol);
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

//----------------------------------------------------------------------------------------------
TPZMaterial *Hexaedro(int grau){
  //Problema teste do Cedric <=> problema teste no paper A. Coutinho e paper Zhang, Yu, Chang levado para 3D
  // e teste no paper de Peyrard and Villedieu
  //  CriacaoDeNos(8,hexaedro);
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
  int index;
  TPZGeoEl *elgc3d = gmesh->CreateGeoElement(ECube,nodes,1,index,1);
  //construtor descont�nuo

  //int interfdim = 2;
  //TPZCompElDisc::gInterfaceDimension = interfdim;
  gmesh->BuildConnectivity();
  int nummat = 1;
	TPZArtDiffType artdiff = LeastSquares_AD;
 // char *artdiff = "LS";
  int nivel;
  cout << "\nmain::Divisao Nivel final da malha ? : ";
  cin >> nivel;
  REAL cfl,delta_x,delta_t,delta,gama;//,maxflinha;

  cfl = ( 1./(2.0*(REAL)grau+1.0) );
  delta_x =  ( 1.0 / pow((REAL)2.0,(REAL)nivel) );
  delta_t = cfl*delta_x;//delta_t � <= que este valor
  //calculando novos valores
  delta_t = delta_x*cfl;
  delta =  (10./3.)*cfl*cfl - (2./3.)*cfl + 1./10.;
  gama = 1.4;
  cout << "\nDominio [0,1]x[0,1]"
       << "\nMax df/dx (desconhecido) = 1.0"
       << "\nCFL = " << cfl
       << "\ndelta otimo = " << delta
       << "\nDelta x = " << delta_x
       << "\ndelta t = " << delta_t
       << "\ndiffusao = " << artdiff << endl;

  int dim = 3;
  TPZMaterial *mat = new TPZEulerConsLaw(nummat,delta_t,gama,dim,artdiff);
  TPZAutoPointer<TPZMaterial> matauto(mat);

  //((TPZConservationLaw *)mat)->SetIntegDegree(grau);
  DebugStop(); //mat->SetForcingFunction(Function);
  cmesh->InsertMaterialObject(matauto);

  //condi��es de contorno
  TPZBndCond *bc;
  //REAL big = 1.e12;
  TPZFMatrix<REAL> val1(5,5,0.),val2(5,1,0.);

  //CC FACE 20: parede
  val1.Zero();
  val2.Zero();
  TPZGeoElBC(elgc3d,20,-1);
  bc = matauto->CreateBC(matauto,-1,5,val1,val2);
  cmesh->InsertMaterialObject(bc);

  //CC FACE 21: Dirichlet
  val1.Zero();
  val2.Zero();
  REAL ro = 1.7;
  REAL u = 2.61934;
  REAL v = 0.50632;
  REAL w = 0.0;
  REAL p = 1.52819;
  REAL vel2 = u*u + v*v + w*w;
  val2(0,0) = ro;
  val2(1,0) = ro * u;
  val2(2,0) = ro * v;
  val2(3,0) = ro * w;
  val2(4,0) = p/(gama-1.0) + 0.5 * ro * vel2;
  TPZGeoElBC(elgc3d,21,-2);
  bc = matauto->CreateBC(matauto,-2,3,val1,val2);
  cmesh->InsertMaterialObject(bc);

  //CC FACE 22  DIREITA : livre
  val1.Zero();
  val2.Zero();
  TPZGeoElBC(elgc3d,22,-3);
  bc = matauto->CreateBC(matauto,-3,4,val1,val2);
  cmesh->InsertMaterialObject(bc);

  //CC FACE 23: POSTERIOR : PAREDE = wall
  val1.Zero();
  val2.Zero();
  TPZGeoElBC(elgc3d,23,-4);
  bc = matauto->CreateBC(matauto,-4,5,val1,val2);//CC MISTA
  cmesh->InsertMaterialObject(bc);

  //CC FACE ESQUERDA 24: Dirichlet
  val1.Zero();
  val2.Zero();
  ro = 1.0;
  u = 2.9;
  v = 0.0;
  w = 0.0;
  p = 0.714286;
  vel2 = u*u+v*v+w*w;
  val2(0,0) = ro;
  val2(1,0) = ro * u;
  val2(2,0) = ro * v;
  val2(3,0) = ro * w;
  val2(4,0) = p/(gama-1.0) +  0.5 * ro * vel2;
  TPZGeoElBC(elgc3d,24,-5);
  bc = matauto->CreateBC(matauto,-5,3,val1,val2);
  cmesh->InsertMaterialObject(bc);

  //CC FACE 25 SUPERIOR : PAREDE
  val1.Zero();
  val2.Zero();
  TPZGeoElBC(elgc3d,25,-6);
  bc = matauto->CreateBC(matauto,-6,5,val1,val2);//CC MISTA
  cmesh->InsertMaterialObject(bc);

  cout << endl;
  cmesh->AutoBuild();
  return mat;
}

//----------------------------------------------------------------------------------------------
TPZMaterial *ProblemaT2D(int grau){
  //Teste no paper de A. Coutinho e primeiro problema teste na tese de Jorge Calle
  //teste do papern Zhang, Yu, Chang e teste no paper de Peyrard and Villedieu
  //  CriacaoDeNos(4,quadrilatero);//para formar dois triangulos
  //elemento de volume
  TPZVec<int> nodes;
  nodes.Resize(3);
  nodes[0] = 0;
  nodes[1] = 2;
  nodes[2] = 3;
  int index;
  TPZGeoEl *elgt2d0 = gmesh->CreateGeoElement(ETriangle,nodes,1,index,1);
  nodes[0] = 0;
  nodes[1] = 1;
  nodes[2] = 2;
  TPZGeoEl *elgt2d1 = gmesh->CreateGeoElement(ETriangle,nodes,1,index,1);

//  int interfdim = 1;
//  TPZCompElDisc::gInterfaceDimension = interfdim;
  gmesh->BuildConnectivity();
  int nummat = 1;
	TPZArtDiffType artdiff = LeastSquares_AD;
  int nivel;
  cout << "\nmain::Divisao Nivel final da malha ? : ";
  cin >> nivel;
  REAL cfl = ( 1./(2.0*(REAL)grau+1.0) );///0.5;
  REAL delta_x =  ( 1.0 / pow((REAL)2.0,(REAL)nivel) );//0.5;
  REAL delta_t = cfl*delta_x;//delta_t � <= que este valor
  //calculando novos valores
  delta_t = delta_x*cfl;
  REAL delta =  (10./3.)*cfl*cfl - (2./3.)*cfl + 1./10.;
  gama = 1.4;
  cout << "\nDominio [0,1]x[0,1]"
       << "\nMax df/dx (desconhecido) = 1.0"
       << "\nCFL = " << cfl
       << "\ndelta otimo = " << delta
       << "\nDelta x = " << delta_x
       << "\ndelta t = " << delta_t
       << "\ndiffusao = " << artdiff << endl;

  int dim = 2;
  TPZMaterial *mat = (TPZEulerConsLaw *) new TPZEulerConsLaw(nummat,delta_t,gama,dim,artdiff);

  DebugStop(); //mat->SetForcingFunction(Function);
  TPZAutoPointer<TPZMaterial> matauto(mat);
  cmesh->InsertMaterialObject(matauto);
  //condi��es de contorno
  TPZBndCond *bc;
  //REAL big = 1.e12;
  TPZFMatrix<REAL> val1(4,4,0.),val2(4,1,0.);

  //TYPE CC
  //0: Dirichlet
  //1: Neumann
  //2: Mista
  //3: Dirichlet -> c�lculo do fluxo
  //4: livre
  //5: parede

  //CC ARESTA INFERIOR
  val1.Zero();
  val2.Zero();
  REAL ro = 1.7;
  REAL u = 2.61934;
  REAL v = 0.50632;
  REAL p = 1.52819;
  REAL vel2 = u*u + v*v;
  val2(0,0) = ro;
  val2(1,0) = ro * u;
  val2(2,0) = ro * v;
  val2(3,0) = p/(gama-1.0) + 0.5 * ro * vel2;
  TPZGeoElBC(elgt2d1,3,-1);
  bc = matauto->CreateBC(matauto,-1,3,val1,val2);
  cmesh->InsertMaterialObject(bc);//bc->SetForcingFunction(Function);

  //CC ARESTA DIREITA : LIVRE
  val1.Zero();
  val2.Zero();
  TPZGeoElBC(elgt2d1,4,-2);
  bc = matauto->CreateBC(matauto,-2,4,val1,val2);
  cmesh->InsertMaterialObject(bc);//bc->SetForcingFunction(Function);

  //CC ARESTA SUPERIOR : PAREDE - WALL
  val1.Zero();
  val2.Zero();
  TPZGeoElBC(elgt2d0,4,-3);
  bc = matauto->CreateBC(matauto,-3,5,val1,val2);
  cmesh->InsertMaterialObject(bc);//bc->SetForcingFunction(Function);

  //CC ARESTA ESQUERDA
  val1.Zero();
  val2.Zero();
  ro = 1.0;
  u = 2.9;
  v = 0.0;
  p = 0.714286;
  val2(0,0) = ro;
  val2(1,0) = ro * u;
  val2(2,0) = ro * v;
  vel2 = u*u+v*v;
  val2(3,0) = p/(gama-1.0) +  0.5 * ro * vel2;
  TPZGeoElBC(elgt2d0,5,-4);
  bc = matauto->CreateBC(matauto,-4,3,val1,val2);
  cmesh->InsertMaterialObject(bc);//bc->SetForcingFunction(Function);

  cout << endl;
  cmesh->AutoBuild();

  return mat;
}

using namespace pzgeom;
using namespace pzrefine;
//----------------------------------------------------------------------------------------------
TPZMaterial *ProblemaQ2D1El(int grau){
  //Teste no paper de A. Coutinho e primeiro problema teste na tese de Jorge Calle
  //teste do papern Zhang, Yu, Chang  e teste no paper de Peyrard and Villedieu
  CriacaoDeNos(4,quadrado);
  //elemento de volume
  TPZVec<int> nodes;
  int index;
  nodes.Resize(4);
  nodes[0] = 0;
  nodes[1] = 1;
  nodes[2] = 2;
  nodes[3] = 3;
  TPZGeoEl *elgq2d = gmesh->CreateGeoElement(EQuadrilateral,nodes,1,index);

  //construtor descont�nuo
	cmesh->SetAllCreateFunctionsDiscontinuous();

//  int interfdim = 1;
//  TPZCompElDisc::gInterfaceDimension = interfdim;
  gmesh->BuildConnectivity();
  int nummat = 1;
	TPZArtDiffType artdiff = LeastSquares_AD;
  //  int nivel = 10;
  //cout << "\nmain::Divisao Nivel final da malha ? : ";
  //cin >> nivel;
  REAL cfl = ( 1./(2.0*(REAL)grau+1.0) );
  REAL delta_x = 1.0;//( 1.0 / pow((REAL)2.0,(REAL)nivel ) );
  REAL delta_t = cfl*delta_x;//delta_t � <= que este valor
  //calculando novos valores
  delta_t = delta_x*cfl;
  REAL delta =  (10./3.)*cfl*cfl - (2./3.)*cfl + 1./10.;
  gama = 1.4;
  cout << "\nDominio [0,1]x[0,1]"
       << "\nMax df/dx (desconhecido) = 1.0"
       << "\nCFL = " << cfl
       << "\ndelta otimo = " << delta
       << "\nDelta x = " << delta_x
       << "\ndelta t = " << delta_t
       << "\ndiffusao = " << artdiff
       << "\ndelta aproximado = " << delta << endl;

  int dim = 2;
  TPZMaterial *mat = new TPZEulerConsLaw(nummat,delta_t,gama,dim,artdiff);
  DebugStop(); //mat->SetForcingFunction(Function);
  TPZAutoPointer<TPZMaterial> matauto(mat);
  cmesh->InsertMaterialObject(matauto);

  //condi��es de contorno
  TPZBndCond *bc;
  //REAL big = 1.e12;
  TPZFMatrix<REAL> val1(4,4),val2(4,1);

  //CC ARESTA INFERIOR
  val1.Zero();
  val2.Zero();
  REAL ro = 1.7;
  REAL u = 2.61934;
  REAL v = 0.50632;
  REAL p = 1.52819;
  REAL vel2 = u*u + v*v;
  val2(0,0) = ro;
  val2(1,0) = ro * u;
  val2(2,0) = ro * v;
  val2(3,0) = p/(gama-1.0) + 0.5 * ro * vel2;
  TPZGeoElBC(elgq2d,4,-1);
  bc = matauto->CreateBC(matauto,-1,3,val1,val2);//bc->SetForcingFunction(Function);
  cmesh->InsertMaterialObject(bc);

  //CC ARESTA DIREITA : LIVRE
  val1.Zero();
  val2.Zero();
  ro = 1.0;
  u = 2.9;
  v = 0.0;
  p = 0.714286;
  vel2 = u*u+v*v;
  val2(0,0) = ro;
  val2(1,0) = ro * u;
  val2(2,0) = ro * v;
  val2(3,0) = p/(gama-1.0) +  0.5 * ro * vel2;
  TPZGeoElBC(elgq2d,5,-2);
  bc = matauto->CreateBC(matauto,-2,4,val1,val2);//bc->SetForcingFunction(Function);
  cmesh->InsertMaterialObject(bc);

  //CC ARESTA SUPERIOR : PAREDE - wall
  val1.Zero();
  val2.Zero();
  TPZGeoElBC(elgq2d,6,-3);
  bc = matauto->CreateBC(matauto,-3,5,val1,val2);//bc->SetForcingFunction(Function);
  cmesh->InsertMaterialObject(bc);

  //CC ARESTA ESQUERDA
  val1.Zero();
  val2.Zero();
  ro = 1.0;
  u = 2.9;
  v = 0.0;
  p = 0.714286;
  val2(0,0) = ro;
  val2(1,0) = ro * u;
  val2(2,0) = ro * v;
  vel2 = u*u+v*v;
  val2(3,0) = p/(gama-1.0) +  0.5 * ro * vel2;
  TPZGeoElBC(elgq2d,7,-4);
  bc = matauto->CreateBC(matauto,-4,3,val1,val2);//bc->SetForcingFunction(Function);
  cmesh->InsertMaterialObject(bc);

  cout << endl;
  cmesh->AutoBuild();

  return mat;
}

//----------------------------------------------------------------------------------------------
TPZMaterial *TresTriangulos(int grau){
  //Teste no paper de A. Coutinho e primeiro problema teste na tese de Jorge Calle
  //teste do papern Zhang, Yu, Chang e teste no paper de Peyrard and Villedieu
  CriacaoDeNos(5,quadrilatero2);//para formar dois triangulos
  //elemento de volume
  TPZVec<int> nodes;
  int index;
  nodes.Resize(3);
  nodes[0] = 0;
  nodes[1] = 1;
  nodes[2] = 3;
  TPZGeoEl *elgt2d0 = gmesh->CreateGeoElement(ETriangle,nodes,1,index);
  nodes[0] = 1;
  nodes[1] = 2;
  nodes[2] = 4;
  TPZGeoEl *elgt2d1 = gmesh->CreateGeoElement(ETriangle,nodes,1,index);
  nodes[0] = 1;
  nodes[1] = 4;
  nodes[2] = 3;
  TPZGeoEl *elgt2d2 = gmesh->CreateGeoElement(ETriangle,nodes,1,index);
	cmesh->SetAllCreateFunctionsDiscontinuous();
//  int interfdim = 1;
//  TPZCompElDisc::gInterfaceDimension = interfdim;
  gmesh->BuildConnectivity();
  int nummat = 1;
//   int nivel;

	TPZArtDiffType artdiff = LeastSquares_AD;
	//cout << "\nmain::Divisao Nivel final da malha ? : ";
  //cin >> nivel;
//   REAL cfl = ( 1./(2.0*(REAL)grau+1.0) );//0.5;
//   REAL delta_x =  ( 1.0 / pow((REAL)2.0,(REAL)nivel) );//0.5;
//   REAL delta_t = cfl*delta_x;//delta_t � <= que este valor
  //calculando novos valores
//   delta_t = delta_x*cfl;
//   REAL delta =  (10./3.)*cfl*cfl - (2./3.)*cfl + 1./10.;
  gama = 1.4;
//   cout << "\nDominio [0,1]x[0,1]"
//        << "\nMax df/dx (desconhecido) = 1.0"
//        << "\nCFL = " << cfl
//        << "\ndelta otimo = " << delta
//        << "\nDelta x = " << delta_x
//        << "\ndelta t = " << delta_t
//        << "\ndiffusao = " << artdiff << endl;

  REAL delta_t = 0.0;
  int dim = 2;
  TPZMaterial *mat = (TPZEulerConsLaw *) new TPZEulerConsLaw(nummat,delta_t,gama,dim,artdiff);

  DebugStop(); //mat->SetForcingFunction(Function);
  TPZAutoPointer<TPZMaterial> matauto(mat);
  cmesh->InsertMaterialObject(matauto);
  //condi��es de contorno
  TPZBndCond *bc;
  TPZFMatrix<REAL> val1(4,4,0.),val2(4,1,0.);

  //CC ARESTA INFERIOR: PAREDE
  val1.Zero();
  val2.Zero();
  TPZGeoElBC(elgt2d0,3,-1);
  TPZGeoElBC(elgt2d1,3,-1);
  bc = matauto->CreateBC(matauto,-1,5,val1,val2);
  cmesh->InsertMaterialObject(bc);

  //CC ARESTA DIREITA : LIVRE
  val1.Zero();
  val2.Zero();
  TPZGeoElBC(elgt2d1,4,-2);
  bc = matauto->CreateBC(matauto,-2,4,val1,val2);
  cmesh->InsertMaterialObject(bc);

  //CC ARESTA SUPERIOR :
  val1.Zero();
  val2.Zero();
  REAL ro = 1.7;
  REAL u = 2.61934;
  REAL v = -0.50632;
  REAL p = 1.52819;
  REAL vel2 = u*u + v*v;
  val2(0,0) = ro;
  val2(1,0) = ro * u;
  val2(2,0) = ro * v;
  val2(3,0) = p/(gama-1.0) + 0.5 * ro * vel2;
  TPZGeoElBC(elgt2d2,4,-3);
  bc = matauto->CreateBC(matauto,-3,3,val1,val2);
  cmesh->InsertMaterialObject(bc);

  //CC ARESTA ESQUERDA
  val1.Zero();
  val2.Zero();
  ro = 1.0;
  u = 2.9;
  v = 0.0;
  p = 0.714286;
  val2(0,0) = ro;
  val2(1,0) = ro * u;
  val2(2,0) = ro * v;
  vel2 = u*u+v*v;
  val2(3,0) = p/(gama-1.0) +  0.5 * ro * vel2;
  TPZGeoElBC(elgt2d0,5,-4);
  bc = matauto->CreateBC(matauto,-4,3,val1,val2);
  cmesh->InsertMaterialObject(bc);

  cout << endl;
  cmesh->AutoBuild();

  return mat;
}

//----------------------------------------------------------------------------------------------
TPZMaterial *TresPrismas(int grau){
  //Problema teste do Cedric <=> problema teste no paper A. Coutinho e paper Zhang, Yu, Chang levado para 3D
  // e teste no paper de Peyrard and Villedieu
  //  CriacaoDeNos(10,tresprismas);
  //elemento de volume
  TPZVec<int> nodes;
  nodes.Resize(6);
  nodes[0] = 0;
  nodes[1] = 1;
  nodes[2] = 3;
  nodes[3] = 5;
  nodes[4] = 6;
  nodes[5] = 8;
  int index;
  TPZGeoEl *elg1 = gmesh->CreateGeoElement(EPrisma,nodes,1,index);
  nodes[0] = 1;
  nodes[1] = 2;
  nodes[2] = 4;
  nodes[3] = 6;
  nodes[4] = 7;
  nodes[5] = 9;
  TPZGeoEl *elg2 = gmesh->CreateGeoElement(EPrisma,nodes,1,index);
  nodes[0] = 3;
  nodes[1] = 1;
  nodes[2] = 4;
  nodes[3] = 8;
  nodes[4] = 6;
  nodes[5] = 9;
  TPZGeoEl *elg3 = gmesh->CreateGeoElement(EPrisma,nodes,1,index);

//  int interfdim = 2;
//  TPZCompElDisc::gInterfaceDimension = interfdim;
  gmesh->BuildConnectivity();
  int nummat = 1,nivel;
	TPZArtDiffType artdiff = LeastSquares_AD;
	cout << "\nmain::Divisao Nivel final da malha ? : ";
  cin >> nivel;
  REAL cfl,delta_x,delta_t,delta,gama;//,maxflinha;

  cfl = ( 1./(2.0*(REAL)grau+1.0) );
  delta_x =  ( 1.0 / pow((REAL)2.0,(REAL)nivel) );
  delta_t = cfl*delta_x;//delta_t � <= que este valor
  delta =  (10./3.)*cfl*cfl - (2./3.)*cfl + 1./10.;
  gama = 1.4;
  cout << "\nDominio [0,1]x[0,1]"
       << "\nMax df/dx (desconhecido) = 1.0"
       << "\nCFL = " << cfl
       << "\ndelta otimo = " << delta
       << "\nDelta x = " << delta_x
       << "\ndelta t = " << delta_t
       << "\ndiffusao = " << artdiff << endl;

  int dim = 3;
  TPZMaterial *mat = (TPZEulerConsLaw *) new TPZEulerConsLaw(nummat,delta_t,gama,dim,artdiff);

  //((TPZConservationLaw *)mat)->SetIntegDegree(grau);
  DebugStop(); //mat->SetForcingFunction(Function);
  TPZAutoPointer<TPZMaterial> matauto(mat);
  cmesh->InsertMaterialObject(matauto);

  //condi��es de contorno
  TPZBndCond *bc;
  TPZFMatrix<REAL> val1(5,5,0.),val2(5,1,0.);

  //CC parede
  val1.Zero();
  val2.Zero();
  TPZGeoElBC(elg1,15,-1);
  TPZGeoElBC(elg2,15,-1);
  TPZGeoElBC(elg3,15,-1);
  TPZGeoElBC(elg1,16,-1);
  TPZGeoElBC(elg2,16,-1);
  TPZGeoElBC(elg1,19,-1);
  TPZGeoElBC(elg2,19,-1);
  TPZGeoElBC(elg3,19,-1);
  bc = matauto->CreateBC(matauto,-1,5,val1,val2);
  cmesh->InsertMaterialObject(bc);

  //CC Dirichlet
  val1.Zero();
  val2.Zero();
  REAL ro = 1.7;
  REAL u = 2.61934;
  REAL v = -0.50632;
  REAL w = 0.0;
  REAL p = 1.52819;
  REAL vel2 = u*u + v*v + w*w;
  val2(0,0) = ro;
  val2(1,0) = ro * u;
  val2(2,0) = ro * v;
  val2(3,0) = ro * w;
  val2(4,0) = p/(gama-1.0) + 0.5 * ro * vel2;
  TPZGeoElBC(elg3,18,-2);
  bc = matauto->CreateBC(matauto,-2,3,val1,val2);
  cmesh->InsertMaterialObject(bc);

  //CC Dirichlet
  val1.Zero();
  val2.Zero();
  ro = 1.0;
  u = 2.9;
  v = 0.0;
  w = 0.0;
  p = 0.714286;
  vel2 = u*u+v*v+w*w;
  val2(0,0) = ro;
  val2(1,0) = ro * u;
  val2(2,0) = ro * v;
  val2(3,0) = ro * w;
  val2(4,0) = p/(gama-1.0) +  0.5 * ro * vel2;
  TPZGeoElBC(elg1,18,-3);
  bc = matauto->CreateBC(matauto,-3,3,val1,val2);
  cmesh->InsertMaterialObject(bc);

  //CC livre
  val1.Zero();
  val2.Zero();
  TPZGeoElBC(elg2,17,-4);
  bc = matauto->CreateBC(matauto,-4,4,val1,val2);
  cmesh->InsertMaterialObject(bc);

  cout << endl;
  cmesh->AutoBuild();
  return mat;
}

//----------------------------------------------------------------------------------------------
TPZMaterial *FluxConst3D(int grau){
  //Problema teste do Cedric <=> problema teste no paper A. Coutinho e paper Zhang, Yu, Chang levado para 3D
  // e teste no paper de Peyrard and Villedieu
  //  CriacaoDeNos(8,hexaedro1);
  //elemento de volume
  TPZVec<int> nodes;
  int index;
  nodes.Resize(8);
  nodes[0] = 0;
  nodes[1] = 1;
  nodes[2] = 2;
  nodes[3] = 3;
  nodes[4] = 4;
  nodes[5] = 5;
  nodes[6] = 6;
  nodes[7] = 7;
  TPZGeoEl *elgc3d = gmesh->CreateGeoElement(ECube,nodes,1,index);
  //construtor descont�nuo
	cmesh->SetAllCreateFunctionsDiscontinuous();
//  int interfdim = 2;
//  TPZCompElDisc::gInterfaceDimension = interfdim;
  gmesh->BuildConnectivity();
  int nummat = 1,nivel;
	TPZArtDiffType artdiff = LeastSquares_AD;
  cout << "\nmain::Divisao Nivel final da malha ? : ";
  cin >> nivel;
  REAL cfl,delta_x,delta_t,delta,gama;//,maxflinha;

  cfl = ( 1./(2.0*(REAL)grau+1.0) );
  delta_x =  ( 1.0 / pow((REAL)2.0,(REAL)nivel) );
  delta_t = cfl*delta_x;//delta_t � <= que este valor
  //calculando novos valores
  delta_t = delta_x*cfl;
  delta =  (10./3.)*cfl*cfl - (2./3.)*cfl + 1./10.;
  gama = 1.4;
  cout << "\nDominio [0,1]x[0,1]"
       << "\nMax df/dx (desconhecido) = 1.0"
       << "\nCFL = " << cfl
       << "\ndelta otimo = " << delta
       << "\nDelta x = " << delta_x
       << "\ndelta t = " << delta_t
       << "\ndiffusao = " << artdiff << endl;

  int dim = 3;
  TPZMaterial *mat = (TPZEulerConsLaw *) new TPZEulerConsLaw(nummat,delta_t,gama,dim,artdiff);

  //((TPZConservationLaw *)mat)->SetIntegDegree(grau);
  DebugStop(); //mat->SetForcingFunction(Function);
  TPZAutoPointer<TPZMaterial> matauto(mat);
  cmesh->InsertMaterialObject(matauto);

  //condi��es de contorno
  TPZBndCond *bc;
  TPZFMatrix<REAL> val1(5,5,0.),val2(5,1,0.);

  //CC FACE: parede
  val1.Zero();
  val2.Zero();
  TPZGeoElBC(elgc3d,20,-1);
  TPZGeoElBC(elgc3d,21,-1);
  TPZGeoElBC(elgc3d,23,-1);
  TPZGeoElBC(elgc3d,25,-1);
  bc = matauto->CreateBC(matauto,-1,5,val1,val2);
  cmesh->InsertMaterialObject(bc);

  //CC FACE : Dirichlet
  val1.Zero();
  val2.Zero();
  REAL ro = 1.0;
  REAL u = 2.9;
  REAL v = 0.0;
  REAL w = 0.0;
  REAL p = 0.714286;
  REAL vel2 = u*u+v*v+w*w;
  val2(0,0) = ro;
  val2(1,0) = ro * u;
  val2(2,0) = ro * v;
  val2(3,0) = ro * w;
  val2(4,0) = p/(gama-1.0) +  0.5 * ro * vel2;
  TPZGeoElBC(elgc3d,24,-2);
  //TPZGeoElBC(elgc3d,22,-2,*gmesh);
  bc = matauto->CreateBC(matauto,-2,3,val1,val2);
  cmesh->InsertMaterialObject(bc);

  //CC FACE 22  DIREITA : livre
  val1.Zero();
  val2.Zero();
  TPZGeoElBC(elgc3d,22,-3);
  bc = matauto->CreateBC(matauto,-3,4,val1,val2);
  cmesh->InsertMaterialObject(bc);

  cout << endl;
  cmesh->AutoBuild();
  return mat;
}


//----------------------------------------------------------------------------------------------
TPZMaterial *FluxConst2D(int grau){

  //  CriacaoDeNos(4,quadrilatero);
  //elemento de volume
  TPZVec<int> nodes;
  nodes.Resize(4);
  nodes[0] = 0;
  nodes[1] = 1;
  nodes[2] = 2;
  nodes[3] = 3;
  int index;
  TPZGeoEl *elgq2d = gmesh->CreateGeoElement(EQuadrilateral,nodes,1,index);

//  int interfdim = 1;
//  TPZCompElDisc::gInterfaceDimension = interfdim;
  gmesh->BuildConnectivity();
  int nummat = 1,nivel;
	TPZArtDiffType artdiff = LeastSquares_AD;
  cout << "\nmain::Divisao Nivel final da malha ? : ";
  cin >> nivel;
  //nivel = 2;
  REAL cfl = ( 1./(2.0*(REAL)grau+1.0) );///0.5;
  REAL delta_x =  ( 1.0 / pow((REAL)2.0,(REAL)nivel) );//0.5;
  REAL delta_t = cfl*delta_x;//delta_t � <= que este valor
  REAL delta =  (10./3.)*cfl*cfl - (2./3.)*cfl + 1./10.;
  gama = 1.4;
  cout << "\nDominio [0,1]x[0,1]"
       << "\nMax df/dx (desconhecido) = 1.0"
       << "\nCFL = " << cfl
       << "\ndelta otimo = " << delta
       << "\nDelta x = " << delta_x
       << "\ndelta t = " << delta_t
       << "\ndiffusao = " << artdiff
       << "\ndelta aproximado = " << delta << endl;

  int dim = 2;
  TPZMaterial *mat = (TPZEulerConsLaw *) new TPZEulerConsLaw(nummat,delta_t,gama,dim,artdiff);
  DebugStop(); //mat->SetForcingFunction(Function);
  TPZAutoPointer<TPZMaterial> matauto(mat);
  cmesh->InsertMaterialObject(matauto);

  //condi��es de contorno
  TPZBndCond *bc;
  TPZFMatrix<REAL> val1(4,4),val2(4,1);
  REAL ro,u,v,vel2,p;

  //CC ARESTAS INFERIOR E SUPERIOR
  val1.Zero();
  val2.Zero();
  TPZGeoElBC(elgq2d,4,-1);
  TPZGeoElBC(elgq2d,6,-1);
  bc = matauto->CreateBC(matauto,-1,5,val1,val2);
  cmesh->InsertMaterialObject(bc);

  //CC ARESTA ESQUERDA
  val1.Zero();
  val2.Zero();
  ro = 1.0;
  u = 2.9;
  v = 0.0;
  p = 0.714286;
  vel2 = u*u+v*v;
  val2(0,0) = ro;
  val2(1,0) = ro * u;
  val2(2,0) = ro * v;
  val2(3,0) = p/(gama-1.0) +  0.5 * ro * vel2;
  TPZGeoElBC(elgq2d,7,-2);
  bc = matauto->CreateBC(matauto,-2,3,val1,val2);
  cmesh->InsertMaterialObject(bc);

  //CC ARESTA DIREITA
  val1.Zero();
  val2.Zero();
  TPZGeoElBC(elgq2d,5,-3);
  bc = matauto->CreateBC(matauto,-3,4,val1,val2);
  cmesh->InsertMaterialObject(bc);

  cout << endl;
  cmesh->AutoBuild();

  return mat;
}

//----------------------------------------------------------------------------------------------
TPZMaterial *NoveQuadrilateros(int grau){
  //Teste no paper de A. Coutinho e primeiro problema teste na tese de Jorge Calle
  //teste do papern Zhang, Yu, Chang  e teste no paper de Peyrard and Villedieu
  CriacaoDeNos(15,novequads);
  //elemento de volume
  TPZVec<int> nodes;
  int INCID[9][4] = {{0,1,6,5},{1,2,7,6},{2,3,10,9},{3,4,11,10},{5,6,7,12},{7,2,9,8},{10,11,14,9},{7,8,13,12},{8,9,14,13}};
  nodes.Resize(4);
  TPZVec<TPZGeoEl *> elem(9);
  elem.Resize(9);
  int i,index;
  for(i=0;i<9;i++){
    nodes[0] = INCID[i][0];
    nodes[1] = INCID[i][1];
    nodes[2] = INCID[i][2];
    nodes[3] = INCID[i][3];
    elem[i] = gmesh->CreateGeoElement(EQuadrilateral,nodes,1,index);
  }
  //construtor descont�nuo
	cmesh->SetAllCreateFunctionsDiscontinuous();
//  int interfdim = 1;
 // TPZCompElDisc::gInterfaceDimension = interfdim;
  gmesh->BuildConnectivity();
  int nummat = 1;
  TPZArtDiffType artdiff = LeastSquares_AD;
  REAL delta_t = 0.0;//ser� calculado
  gama = 1.4;//ar
  int dim = 2;
  TPZMaterial *mat = (TPZEulerConsLaw *) new TPZEulerConsLaw(nummat,delta_t,gama,dim,artdiff);
  DebugStop(); //mat->SetForcingFunction(Function);
  TPZAutoPointer<TPZMaterial> matauto(mat);
  cmesh->InsertMaterialObject(matauto);

  //condi��es de contorno
  TPZBndCond *bc;
  TPZFMatrix<REAL> val1(4,4),val2(4,1);

  //CC ARESTA INFERIOR : PAREDE
  val1.Zero();
  val2.Zero();
  TPZGeoElBC((TPZGeoEl  *)elem[0],4,-1);
  TPZGeoElBC((TPZGeoEl  *)elem[1],4,-1);
  TPZGeoElBC((TPZGeoEl  *)elem[2],4,-1);
  TPZGeoElBC((TPZGeoEl  *)elem[3],4,-1);
  bc = matauto->CreateBC(matauto,-1,5,val1,val2);
  cmesh->InsertMaterialObject(bc);

  //CC ARESTA DIREITA : OUTFLOW
  val1.Zero();
  val2.Zero();
  TPZGeoElBC((TPZGeoEl  *)elem[3],5,-2);
  TPZGeoElBC((TPZGeoEl  *)elem[6],5,-2);
  bc = matauto->CreateBC(matauto,-2,4,val1,val2);
  cmesh->InsertMaterialObject(bc);

  //CC ARESTA SUPERIOR : DIRICHLET
  val1.Zero();
  val2.Zero();
  REAL ro = 1.7;
  REAL u = 2.61934;
  REAL v = -0.50632;
  REAL p = 1.52819;
  REAL vel2 = u*u + v*v;
  val2(0,0) = ro;
  val2(1,0) = ro * u;
  val2(2,0) = ro * v;
  val2(3,0) = p/(gama-1.0) + 0.5 * ro * vel2;
  TPZGeoElBC((TPZGeoEl  *)elem[7],6,-3);
  TPZGeoElBC((TPZGeoEl  *)elem[8],6,-3);
  bc = matauto->CreateBC(matauto,-3,3,val1,val2);
  cmesh->InsertMaterialObject(bc);

  //CC ARESTA ESQUERDA : INFLOW
  val1.Zero();
  val2.Zero();
  ro = 1.0;
  u = 2.9;
  v = 0.0;
  p = 0.714286;
  val2(0,0) = ro;
  val2(1,0) = ro * u;
  val2(2,0) = ro * v;
  vel2 = u*u+v*v;
  val2(3,0) = p/(gama-1.0) +  0.5 * ro * vel2;
  TPZGeoElBC((TPZGeoEl  *)elem[0],7,-4);
  TPZGeoElBC((TPZGeoEl  *)elem[4],7,-4);
  bc = matauto->CreateBC(matauto,-4,3,val1,val2);
  cmesh->InsertMaterialObject(bc);

  cout << endl;
  cmesh->AutoBuild();

  return mat;
}
////////////////////////////////////
void Function(TPZVec<REAL> &x,TPZVec<REAL> &result){

  if(problem == -1){
    result.Resize(4);
    for(int i=0;i<4;i++) result[i] = 1.;
    return;
  }

  if(problem == 6){
    result.Resize(5);
    //Condi��o inicial t =  0
    REAL ro = 1.0;
    REAL u = 2.9;// 1.0 , 2.9
    REAL v = 0.0;
    REAL w = 0.0;
    REAL p = 0.714286;// 0.246306 , 2.9
    REAL vel2 = u*u+v*v+w*w;
    result[0] = ro;
    result[1] = ro * u;
    result[2] = ro * v;
    result[3] = ro * w;
    result[4] = p/(gama-1.0) +  0.5 * ro * vel2;
    return;
  }

  if(problem == 0){
    result.Resize(4);
    //Condi��o inicial t =  0
    REAL ro = 1.0;
    REAL u = 2.9;
    REAL v = 0.0;
    REAL p = 0.714286;
    REAL vel2 = u*u + v*v;
    result[0] = ro;
    result[1] = ro * u;
    result[2] = ro * v;
    result[3] = p/(gama-1.0) + 0.5 * ro * vel2;
    return;
  }

  if(problem == 1){
    result.Resize(5);
    //Condi��o inicial t =  0
    REAL ro = 1.0;
    REAL u = 2.9;
    REAL v = 0.0;
    REAL w = 0.0;
    REAL p = 0.714286;
    REAL vel2 = u*u+v*v+w*w;
    result[0] = ro;
    result[1] = ro * u;
    result[2] = ro * v;
    result[3] = ro * w;
    result[4] = p/(gama-1.0) +  0.5 * ro * vel2;
    return;
  }
  if(problem == 2){
    result.Resize(5);
    //Condi��o inicial t =  0
    REAL ro = 1.0;
    REAL u = 2.9;// 1.0 , 2.9
    REAL v = 0.0;
    REAL w = 0.0;
    REAL p = 0.714286;// 0.246306 , 2.9
    REAL vel2 = u*u+v*v+w*w;
    result[0] = ro;
    result[1] = ro * u;
    result[2] = ro * v;
    result[3] = ro * w;
    result[4] = p/(gama-1.0) +  0.5 * ro * vel2;
    return;

    if(0){
      REAL ro = 1.0;
      REAL u = 2.9;
      REAL v = 0.0;
      REAL w = 0.0;
      REAL p = 0.714286;
      REAL vel2 = u*u+v*v+w*w;
      result[0] = ro;
      result[1] = ro * u;
      result[2] = ro * v;
      result[3] = ro * w;
      result[4] = p/(gama-1.0) +  0.5 * ro * vel2;
      return;
    }
  }
  if(problem == 3 || problem == 7){
    result.Resize(4);
    //Condi��o inicial t =  0
    REAL ro = 1.0;
    REAL u = 2.9;
    REAL v = 0.0;
    REAL p = 0.714286;
    REAL vel2 = u*u+v*v;
    result[0] = ro;
    result[1] = ro * u;
    result[2] = ro * v;
    result[3] = p/(gama-1.0) +  0.5 * ro * vel2;
    return;
  }
  if(problem == 4){
    result.Resize(4);
    //Condi��o inicial t =  0
    REAL ro = 1.0;
    REAL u = 1.0;
    REAL v = 0.0;
    REAL p = 0.246306;
    REAL vel2 = u*u+v*v;
    result[0] = ro;
    result[1] = ro * u;
    result[2] = ro * v;
    result[3] = p/(gama-1.0) +  0.5 * ro * vel2;
    return;
  }
}

//----------------------------------------------------------------------------------------------
TPZMaterial *NoveCubos(int grau){
  //Problema teste do Cedric <=> problema teste no paper A. Coutinho e paper Zhang, Yu, Chang levado para 3D
  // e teste no paper de Peyrard and Villedieu
   CriacaoDeNos(30,novecubos);
  //elemento de volume
  TPZVec<int> nodes;
  int INCID[9][8] = {{0,1,6,5,15,16,21,20},{1,2,7,6,16,17,22,21},{2,3,10,9,17,18,25,24},{3,4,11,10,18,19,26,25},
		     {5,6,7,12,20,21,22,27},{7,2,9,8,22,17,24,23},{10,11,14,9,25,26,29,24},{7,8,13,12,22,23,28,27},
		     {8,9,14,13,23,24,29,28}};
  nodes.Resize(8);
  TPZVec<TPZGeoEl *> elem(9);
  elem.Resize(9);
  int i,index;
  for(i=0;i<9;i++){
    nodes[0] = INCID[i][0];
    nodes[1] = INCID[i][1];
    nodes[2] = INCID[i][2];
    nodes[3] = INCID[i][3];
    nodes[4] = INCID[i][4];
    nodes[5] = INCID[i][5];
    nodes[6] = INCID[i][6];
    nodes[7] = INCID[i][7];
    elem[i] = gmesh->CreateGeoElement(ECube,nodes,1,index);
  }

  //elemento de volume descont�nuo
	cmesh->SetAllCreateFunctionsDiscontinuous();
//  int interfdim = 2;
//  TPZCompElDisc::gInterfaceDimension = interfdim;
  gmesh->BuildConnectivity();
  int nummat = 1;
	TPZArtDiffType artdiff = LeastSquares_AD;
  REAL delta_t = 0.0;//ser� calculado
  gama = 1.4;//ar
  int dim = 3;
  TPZMaterial *mat = (TPZEulerConsLaw *) new TPZEulerConsLaw(nummat,delta_t,gama,dim,artdiff);

  //((TPZConservationLaw *)mat)->SetIntegDegree(grau);
  DebugStop(); //mat->SetForcingFunction(Function);
  TPZAutoPointer<TPZMaterial> matauto(mat);
  cmesh->InsertMaterialObject(matauto);

  //condi��es de contorno
  TPZBndCond *bc;
  TPZFMatrix<REAL> val1(5,5,0.),val2(5,1,0.);

  //CC ARESTA INFERIOR : PAREDE
  val1.Zero();
  val2.Zero();
  TPZGeoElBC((TPZGeoEl  *)elem[0],20,-1);
  TPZGeoElBC((TPZGeoEl  *)elem[1],20,-1);
  TPZGeoElBC((TPZGeoEl  *)elem[2],20,-1);
  TPZGeoElBC((TPZGeoEl  *)elem[3],20,-1);
  TPZGeoElBC((TPZGeoEl  *)elem[4],20,-1);
  TPZGeoElBC((TPZGeoEl  *)elem[5],20,-1);
  TPZGeoElBC((TPZGeoEl  *)elem[6],20,-1);
  TPZGeoElBC((TPZGeoEl  *)elem[7],20,-1);
  TPZGeoElBC((TPZGeoEl  *)elem[8],20,-1);
  TPZGeoElBC((TPZGeoEl  *)elem[0],25,-1);
  TPZGeoElBC((TPZGeoEl  *)elem[1],25,-1);
  TPZGeoElBC((TPZGeoEl  *)elem[2],25,-1);
  TPZGeoElBC((TPZGeoEl  *)elem[3],25,-1);
  TPZGeoElBC((TPZGeoEl  *)elem[4],25,-1);
  TPZGeoElBC((TPZGeoEl  *)elem[5],25,-1);
  TPZGeoElBC((TPZGeoEl  *)elem[6],25,-1);
  TPZGeoElBC((TPZGeoEl  *)elem[7],25,-1);
  TPZGeoElBC((TPZGeoEl  *)elem[8],25,-1);
  TPZGeoElBC((TPZGeoEl  *)elem[0],21,-1);
  TPZGeoElBC((TPZGeoEl  *)elem[1],21,-1);
  TPZGeoElBC((TPZGeoEl  *)elem[2],21,-1);
  TPZGeoElBC((TPZGeoEl  *)elem[3],21,-1);
  bc = matauto->CreateBC(matauto,-1,5,val1,val2);
  cmesh->InsertMaterialObject(bc);

  //CC ARESTA DIREITA : OUTFLOW
  val1.Zero();
  val2.Zero();
  TPZGeoElBC((TPZGeoEl  *)elem[3],22,-2);
  TPZGeoElBC((TPZGeoEl  *)elem[6],22,-2);
  bc = matauto->CreateBC(matauto,-2,4,val1,val2);
  cmesh->InsertMaterialObject(bc);

  //CC ARESTA SUPERIOR : DIRICHLET
  val1.Zero();
  val2.Zero();
  REAL ro = 1.7;
  REAL u = 2.61934;
  REAL v = -0.50632;
  REAL w = 0.0;
  REAL p = 1.52819;
  REAL vel2 = u*u + v*v + w*w;
  val2(0,0) = ro;
  val2(1,0) = ro * u;
  val2(2,0) = ro * v;
  val2(3,0) = ro * w;
  val2(4,0) = p/(gama-1.0) + 0.5 * ro * vel2;
  TPZGeoElBC((TPZGeoEl  *)elem[7],23,-3);
  TPZGeoElBC((TPZGeoEl  *)elem[8],23,-3);
  bc = matauto->CreateBC(matauto,-3,3,val1,val2);
  cmesh->InsertMaterialObject(bc);

  //CC ARESTA ESQUERDA : INFLOW
  val1.Zero();
  val2.Zero();
  ro = 1.0;
  u = 2.9;
  v = 0.0;
  w = 0.0;
  p = 0.714286;
  vel2 = u*u+v*v+w*w;
  val2(0,0) = ro;
  val2(1,0) = ro * u;
  val2(2,0) = ro * v;
  val2(3,0) = ro * w;
  val2(4,0) = p/(gama-1.0) +  0.5 * ro * vel2;
  TPZGeoElBC((TPZGeoEl  *)elem[0],24,-4);
  TPZGeoElBC((TPZGeoEl  *)elem[4],24,-4);
  bc = matauto->CreateBC(matauto,-4,3,val1,val2);
  cmesh->InsertMaterialObject(bc);

  cout << endl;
  cmesh->AutoBuild();
  return mat;
}

void SequenceDivide(int fat[100],int numbel){

  TPZVec<int> csub(0);
  int i,el;
  for(i=0;i<numbel;i++){
    cout << "\nId do elemento geometrico a dividir -> " << fat[i] << endl;
    cmesh->Divide(fat[i],csub,1);
  }
  cout << "Elementos divissiveis:\n";
  numbel = cmesh->ElementVec().NElements();
  for(el=0;el<numbel;el++) {
    TPZCompEl *cpel = cmesh->ElementVec()[el];
    if(cpel && cpel->Type() != EInterface){
      TPZGeoEl *gel = cpel->Reference();
      if(gel) cout << gel->Id() << ",";
    }
  }
}

void SequenceDivide2(){

  TPZVec<int> s(0),s2(0);
  cmesh->Divide(1,s,0);
  int niv = 0;
  if(0){
    while(niv++ < 4){
      cmesh->Divide(s[0],s2,0);
      cmesh->Divide(s[1],s2,0);
      cmesh->Divide(s[3],s2,0);
      cmesh->Divide(s[2],s,0);
    }
    niv = 0;
    cmesh->Divide(2,s,0);
    while(niv++ < 4){
      cmesh->Divide(s[0],s2,0);
      cmesh->Divide(s[3],s2,0);
      cmesh->Divide(s[1],s,0);
    }
  }
  if(0){
    int i;//novecubos
    for(i=0;i<9;i++) cmesh->Divide(i,s,0);
  }
  if(1){//1 cubos
    cmesh->Divide(0,s,0);
    cmesh->Divide(s[7],s,0);//6
    cmesh->Divide(s[5],s,0);//4
    cmesh->Divide(s[5],s,0);//4: essa combina��o da problemas
    //cmesh->Divide(s[4],s,0);
  }
}

void AgrupaList(TPZVec<int> &accumlist,int nivel,int &numaggl){
  //todo elemento deve ser agrupado nem que for para ele mesmo
  cmesh->SetDimModel(meshdim);
  cout << "\n\nmain::AgrupaList para malha 2D\n\n";
  int nel = cmesh->NElements(),i;
  //n�o todo index � sub-elemento
  accumlist.Resize(nel,-1);
  int mdim = cmesh->Dimension();
  int sdim = mdim - 1;
  for(i=0;i<nel;i++){
    TPZCompEl *cel = cmesh->ElementVec()[i];
    if(!cel) continue;
    TPZGeoEl *gel = cel->Reference();
    int type = cel->Type(),eldim = cel->Dimension();
    //agrupando elementos computacionais de volume
    if(type == EInterface) continue;//interface ser� clonada
    if(eldim == sdim) continue;//discontinuo BC ser� clonado
    TPZGeoEl *father = gel->Father();
    if(!father) continue;
    while(father && Nivel(father) != nivel) father = father->Father();//nova
    if (!father) continue;//nova
    //if(Nivel(father) != nivel) continue;//antiga
    int fatid = father->Id();
    accumlist[i] = fatid;
  }
  //reordena a lista por ordem crescente do pai
  TPZVec<int> list(accumlist);
  int j;
  for(i=0;i<nel;i++){
    for(j=i+1;j<nel;j++){
      if(list[i] > list[j]){
	int aux = list[i];
	list[i] = list[j];
	list[j] = aux;
      }
    }
  }
  //conta o n�mero de elementos obtidos por aglomera��o
  numaggl = 0;
  int act2 = -1;
  for(i=0;i<nel;i++){
    int act = list[i];
    if(act == act2) continue;
    for(j=i+1;j<nel;j++){
      int next = list[j];
      if(next != act2){
	numaggl++;
	j = nel;
	act2 = act;
      }
    }
  }
  //reformula o pai de 0 a nmax
  TPZVec<int> newlist(accumlist);
  int newfat = 0;
  for(i=0;i<nel;i++){
    int fatid1 = newlist[i];
    if(fatid1 < 0) continue;
    accumlist[i] = newfat;
    newlist[i] = -1;
    for(j=i+1;j<nel;j++){
      int fatid2 = newlist[j];
      if(fatid2 == fatid1){
	accumlist[j] = newfat;
	newlist[j] = -1;
      }
    }
    newfat++;
  }
  if(newfat != numaggl) cout << "main::AgrupaList n�mero de pais n�o confere\n";
  if(!newfat && !numaggl) cout << "main::AgrupaList lista de elementos aglomerados vacia\n";
}


/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
TPZMaterial *Quadrado(int grau){

  //  CriacaoDeNos(4,quadrado);
  //elemento de volume
  TPZVec<int> nodes;
  nodes.Resize(4);
  TPZGeoEl *elem;
  int index;
  nodes[0] = 0;
  nodes[1] = 1;
  nodes[2] = 2;
  nodes[3] = 3;
  elem = gmesh->CreateGeoElement(EQuadrilateral,nodes,1,index);

  //construtor descont�nuo
	cmesh->SetAllCreateFunctionsDiscontinuous();
//  int interfdim = 1;
//  TPZCompElDisc::gInterfaceDimension = interfdim;
  gmesh->BuildConnectivity();
  //  int nummat = 1;
  //  int dim = 2;
  TPZMaterial *mat;// = new TPZMatHybrid(nummat,dim);
  DebugStop(); //mat->SetForcingFunction(Function);
  TPZAutoPointer<TPZMaterial> matauto(mat);
  cmesh->InsertMaterialObject(matauto);

  //condi��es de contorno
  TPZBndCond *bc;
  TPZFMatrix<REAL> val1(1,1,0.),val2(1,1,0.);

  //CC DE NEUMANN

  TPZGeoElBC(elem,5,-1);
  bc = matauto->CreateBC(matauto,-1,1,val1,val2);
  //  bc->SetForcingFunction(G_Function);
  cmesh->InsertMaterialObject(bc);

  cmesh->AutoBuild();

  return mat;
}

void SetDeltaTime(TPZMaterial *mat,TPZCompMesh *cmesh){

  TPZVec<REAL> x(3,0.0),sol;
  int i,nstate = mat->NStateVariables();
  cout << "main::SetDeltaTime ENTRE PONTO NO DOMINIO\n";
  x[0] = 2.0;
  x[1] = 0.5;
  Function(x,sol);
  REAL prod = 0.0,maxveloc;
  for(i=1;i<nstate-1;i++) prod += sol[i]*sol[i];//(u�+v�+w�)*ro�
  REAL dens2 = sol[0]*sol[0];
  maxveloc = sqrt(prod/dens2);//velocidade
  TPZEulerConsLaw *law = dynamic_cast<TPZEulerConsLaw *>(mat);
  //TPZEulerConsLaw *law = (TPZEulerConsLaw *)(mat);
  REAL press = law->Pressure(sol);
  if(press < 0)
    cout << "main::SetDeltaTime pressao <0 , toma valor absoluto para calculo do som\n";
  REAL sound = sqrt(gama*press/sol[0]);
  maxveloc += sound;
  //REAL deltax = cmesh->DeltaX();
  REAL deltax = cmesh->LesserEdgeOfMesh();
  //REAL deltax = cmesh->MaximumRadiusOfMesh();
  deltaT = CFL*deltax/maxveloc;//static REAL
  cout << "main::SetDeltaTime : " << deltaT << endl;
  law->SetDelta(deltaT);
  int continua;
  cout <<  "\n -> ";
  cin >> continua;

}

//   cout << "\nmain::Divisao Nivel final da malha ? : ";
//   cin >> nivel;
//   REAL cfl,delta_x,delta_t,delta,gama;//,maxflinha;

//   cfl = ( 1./(2.0*(REAL)grau+1.0) );
//   delta_x =  ( 1.0 / pow((REAL)2.0,(REAL)nivel) );
//   delta_t = cfl*delta_x;//delta_t � <= que este valor
//   delta =  (10./3.)*cfl*cfl - (2./3.)*cfl + 1./10.;
//   gama = 1.4;
//   cout << "\nDominio [0,1]x[0,1]"
//        << "\nMax df/dx (desconhecido) = 1.0"
//        << "\nCFL = " << cfl
//        << "\ndelta otimo = " << delta
//        << "\nDelta x = " << delta_x
//        << "\ndelta t = " << delta_t
//        << "\ndiffusao = " << artdiff << endl;


//static TPZRefPattern hexaedro9("hexaedro9subs.in");
//static TPZRefPattern hexaedro8("hexaedro8.in");//hexaedro mestre 8 sub-elementos
//static TPZRefPattern hexaedro2("hexaedro2.in");
//static TPZRefPattern quad4("quadrilatero4.in");//quadrilatero mestre 4 sub-elementos*/
//static TPZRefPattern linha2("linha2.in");//aresta mestre 2 sub-elementos*/

//   if(0){
//     TPZFMatrix<REAL> TS(4,4,0.),D(4,4,0.),TI(4,4),mat(4,4,0.);
//     TS(0,0) = 1.0; TS(0,1) = 1.0; TS(0,2) = 2.0; TS(0,3) = 4.0;
//     TS(1,1) = 1.0; TS(1,2) = 3.0; TS(1,3) = 5.0;
//     TS(2,2) = 1.0; TS(2,3) = 6.0;
//     TS(3,3) = 1.0;
//     TS.Transpose(&TI);
//     D(0,0) = 10.0;
//     D(1,1) = 20.0;
//     D(2,2) = 30.0;
//     D(3,3) = 40.0;
//     ofstream mats("MATRIZES_TESTES.out");
//     mats << "\n\n\t\t\t* * * MATRIZES TESTES * * *\n\n";
//     TI.Print("\n\n\t\t\t* * * MATRIZES INFERIOR * * *\n\n",mats);
//      D.Print("\n\n\t\t\t* * * MATRIZES DIAGONAL * * *\n\n",mats);
//     TS.Print("\n\n\t\t\t* * * MATRIZES SUPERIOR * * *\n\n",mats);
//     mat = TI * (D * TS);
//     mat.Print("\n\n\t\t\t* * * MATRIZES PRODUTO TI*D*TS * * *\n\n",mats);

//     TPZNonLinMultGridAnalysis an(cmesh);
//     TPZFStructMatrix stiff(cmesh);
//     an.SetStructuralMatrix(stiff);
//     an.Solution().Zero();
//     TPZStepSolver solver;
//     solver.SetDirect(ELDLt);//ELU, ECholesky
//     an.SetSolver(solver);
//     an.Solver().SetMatrix(&mat);
//     TPZFMatrix<REAL> sol(4,1,0.),res(4,1,0.),anres(4,1,0.);
//     sol(0,0) = 1.0;
//     sol(1,0) = 1.0;
//     sol(2,0) = 1.0;
//     sol(3,0) = 1.0;
//     an.CalcResidual(sol,anres,res,an,"LDLt");
//     res.Print("\n\n\t\t\t* * * MATRIZ RES�DUO * * *\n\n",cout);
//     mats.flush();
//     mats.close();

//   }
//     //if(cmesh) delete cmesh;
//   //    if(gmesh) delete gmesh;
//   //  return 1;


