

#include "pzfmatrix.h"
//#include "pzskylmat.h"
#include "pzskylstrmatrix.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzstepsolver.h"
#include "pzelasmat.h"
#include "pzbndcond.h"
#include <stdlib.h>
#include <iostream>
#include "pzvec.h"
#include "pzstack.h"
#include "pzcompel.h"
#include "pzgeoel.h"

//#include "pzelgq2d.h"
//#include "pzelgt2d.h"
#include "pzmat2dlin.h"
#include "pzanalysis.h"
#include "pzmetis.h"
#include "tpznodesetcompute.h"
#include "pzysmp.h"
#include "TPZTimer.h"

#include <stdio.h>
#include <time.h>
#include <math.h>

#include "fast.h"
//#include "pzelct2d.h"
//template<class T>
//class TPZVec;
//#define NOTDEBUG

void InicializarMaterial(TPZCompMesh &cmesh);
void LerMalha(char *arquivo,TPZGeoMesh &mesh);
void UmElemento(TPZGeoMesh &malha);
void forcingfunction(TPZVec<REAL> &ponto, TPZVec<REAL> &force);

int main() {
  
  
//  string filename("CursoAltoDes.in");
//  string filename("matriz1089.in");
//  string filename("matriz12.in");
//  string filename("matriz16000.in");
  string filename("matriz260000.in");
  TPZFMatrix rhs,sol;
  
  TPZMultiTimer timer(5);
  cout << "Leitura\n";
  timer.processName(0) = "Leitura da matriz";
  timer.start(0);
  TPZFYsmpMatrix *faster = ReadMatrix(filename,rhs);
  timer.stop(0);
  
  cout << "Multiplicacao Numero de termos " << faster->NumTerms() << "\n";
  TimeMultiply(faster,timer);
  
  cout << "Jacobi\n";
  SolveJacobi(faster, rhs, 1.e-8,timer);
  
  cout << "SSor\n";
  SolveSSOR(faster,rhs,1.e-8,timer);
  
  cout << "CG\n";
  SolveCG(faster,rhs,1.e-8,timer);
  
  cout << timer;
/*
  TPZFMatrix test(*faster);
  
  cputime = SolveJacobi(&test,rhs,1.e-8);
  
  cout << "Tempo para Jacobi "<< cputime << endl;

  cputime = SolveSSOR(&test,rhs,1.e-8);
  
  cout << "Tempo para SSOR "<< cputime << endl;

  cputime = SolveCG(&test,rhs,1.e-8);
  
  cout << "Tempo para CG "<< cputime << endl;
*/
  
//  TPZStepSolver stepsolve(faster);

  return 0;
  
  std::ofstream file("graphtest.txt");
   //malha geometrica
   TPZGeoMesh *malha = new TPZGeoMesh;
   malha->SetName("Malha gerada por Paulo Lyra, nao tem garantia nemhuma");
   

   UmElemento(*malha);
   //LerMalha("quad_st_800_1R_16X.gri",*malha);

   //ordem de interpolacao
   int ord = 1;
   cout << "Entre ordem 1,2,3,4,5 : ";
//   cin >> ord;
   TPZCompEl::gOrder = ord;
   //construção malha computacional
   TPZVec<int> csub(0);
   TPZManVector<TPZGeoEl *> pv(4);
   int n1=1,level=0; 
   cout << "\nDividir ate nivel ? ";
   int resp = 1;
//   cin >> resp;      
   int nelc = malha->ElementVec().NElements();
   int el;  

   TPZGeoEl *cpel;
   for(el=0;el<malha->ElementVec().NElements();el++) {
     cpel = malha->ElementVec()[el];
     if(cpel && cpel->Level() < resp)
		cpel->Divide(pv);

   }
   //analysis
   TPZCompMesh *malhacomp = new TPZCompMesh(malha);
   InicializarMaterial(*malhacomp);
   malhacomp->AutoBuild();
   malhacomp->AdjustBoundaryElements();
   malhacomp->InitializeBlock();
/*
   TPZNodesetCompute nodeset;
   TPZStack<int> elementgraph,elementgraphindex;
   malhacomp->ComputeElGraph(elementgraph,elementgraphindex);
   int nindep = malhacomp->NIndependentConnects();
   malhacomp->ComputeElGraph(elementgraph,elementgraphindex);
   int nel = elementgraphindex.NElements()-1;
   TPZMetis renum(nel,nindep);
   nodeset.Print(file,elementgraphindex,elementgraph);
   renum.ConvertGraph(elementgraph,elementgraphindex,nodeset.Nodegraph(),nodeset.Nodegraphindex());
//   cout << "nodegraphindex " << nodeset.Nodegraphindex() << endl;
//   cout << "nodegraph " << nodeset.Nodegraph() << endl;
   nodeset.AnalyseGraph();
   nodeset.Print(file);
   TPZStack<int> blockgraph,blockgraphindex;
//   nodeset.BuildNodeGraph(blockgraph,blockgraphindex);
//   nodeset.BuildVertexGraph(blockgraph,blockgraphindex);
   nodeset.BuildElementGraph(blockgraph,blockgraphindex);
   
   file << "Vertex graph\n";
   nodeset.Print(file,blockgraphindex,blockgraph);

*/
   TPZAnalysis an(malhacomp);
   TPZSkylineStructMatrix strmat(malhacomp);
   an.SetStructuralMatrix(strmat);
   
   TPZStepSolver step(strmat.Create());
//   step.SetDirect(ECholesky);   
   an.SetSolver(step);
   TPZMatrixSolver *precond = an.BuildPreconditioner(TPZAnalysis::EElement,false);
   step.SetCG(100,*precond,1.e-10,0);
   an.SetSolver(step);
   delete precond;

   malhacomp->SetName("Malha Computacional :  Connects e Elementos");
   // Posprocessamento
   
   an.Run();
   TPZVec<char *> scalnames(1);
   scalnames[0] = "state";

   TPZVec<char *> vecnames(0);

   char plotfile[] =  "termica.dx";
   an.DefineGraphMesh(2, scalnames, vecnames, plotfile);
   an.PostProcess(1);
   //   an.DefineGraphMesh(2, scalnames, vecnames, pltfile);
   //   an.PostProcess(2);
   delete malhacomp;
   delete malha;
   return 0;
}

void UmElemento(TPZGeoMesh &malha) {
	double coordstore[4][3] = {{0.,0.,0.},{1.,0.,0.},{1.,1.,0.},{0.,1.,0.}};
	// criar um objeto tipo malha geometrica

	// criar quatro nos
	int i,j;
	TPZVec<REAL> coord(3,0.);
	for(i=0; i<4; i++) {
		// initializar as coordenadas do no em um vetor
		for (j=0; j<3; j++) coord[j] = coordstore[i][j];

		// identificar um espaço no vetor onde podemos armazenar
		// este vetor
		int nodeindex = malha.NodeVec ().AllocateNewElement ();

		// initializar os dados do nó
		malha.NodeVec ()[i].Initialize (i,coord,malha);
	}

	// criar um elemento

	// initializar os indices dos nós
	TPZVec<int> indices(4);
	for(i=0; i<4; i++) indices[i] = i;

	// O proprio construtor vai inserir o elemento na malha
//	TPZGeoEl *gel = new TPZGeoElQ2d(0,indices,1,malha);
  int index;
	TPZGeoEl *gel = malha.CreateGeoElement(EQuadrilateral,indices,1,index);


	malha.BuildConnectivity ();

	// Associar o lado de um elemento com uma condicao de contorno
	// Este objeto ira inserir-se automaticamente na malha
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


void LerMalha(char *nome, TPZGeoMesh &grid) {
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

	TPZMat2dLin *meumat = new TPZMat2dLin(1);
	TPZFMatrix xk(1,1,1.),xc(1,2,0.),xf(1,1,1.);
	meumat->SetMaterial (xk,xc,xf);
//	meumat->SetForcingFunction(forcingfunction);
	cmesh.InsertMaterialObject(meumat);

	// inserir a condicao de contorno
	TPZFMatrix val1(1,1,0.),val2(1,1,0.);
	TPZMaterial *bnd = meumat->CreateBC (-4,0,val1,val2);
	cmesh.InsertMaterialObject(bnd);
	bnd = meumat->CreateBC (-3,0,val1,val2);
	cmesh.InsertMaterialObject(bnd);
	bnd = meumat->CreateBC (-2,0,val1,val2);
	cmesh.InsertMaterialObject(bnd);
	bnd = meumat->CreateBC (-1,0,val1,val2);
	cmesh.InsertMaterialObject(bnd);


}

