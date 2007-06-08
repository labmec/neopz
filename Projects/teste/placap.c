

#include "pzfmatrix.h"
#include "pzskylmat.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzgnode.h"
#include "pzsolve.h"
#include "pzelg1d.h"
#include "pzmat1dlin.h"
#include "pzelmat.h"
#include "pzelasmat.h"
#include "pzcompel.h"
#include "pzbndcond.h"
#include "pzstack.h"
#include "pztrnsform.h"
#include <stdlib.h>
#include <iostream.h>
#include "pzvec.h"

#include "pzelgq2d.h"
#include "pzelgt2d.h"
#include "pzmat2dlin.h"
#include "pzanalysis.h"
#include "pzmetis.h"
#include "pzplaca.h"

#include <stdio.h>
#include <time.h>
#include "pzelct2d.h"
template<class T>
class TPZVec;
#define NOTDEBUG

void ExactSol(TPZVec<REAL> &loc,TPZVec<REAL> &val,TPZMatrix &deriv);
TPZMaterial *LerMaterial(char *filename);
int main() {

   //malha geometrica
   TPZGeoMesh *firstmesh = new TPZGeoMesh;
   firstmesh->SetName("Malha Geometrica : Nós e Elementos");
   firstmesh->NodeVec().Resize(4);
   TPZVec<REAL> coord(2);
   coord[0] = 0.;
   coord[1] = 0.;
   //nos geometricos
   firstmesh->NodeVec()[0].Initialize(coord,*firstmesh);
   coord[0] = 1.0;
   firstmesh->NodeVec()[1].Initialize(coord,*firstmesh);
   coord[0] = 1.0;
   coord[1] = 1.0;
   firstmesh->NodeVec()[2].Initialize(coord,*firstmesh);
   coord[0] = 0.0;
   firstmesh->NodeVec()[3].Initialize(coord,*firstmesh);
   TPZVec<int> nodeindexes(4);
   nodeindexes[0] = 0;
   nodeindexes[1] = 1;
   nodeindexes[2] = 2;
   nodeindexes[3] = 3;
    //elementos geometricos
   TPZGeoElQ2d *elg0 = new TPZGeoElQ2d(nodeindexes,1,*firstmesh);
   //Arquivos de saida
	ofstream outgm1("outgm1.dat");
   ofstream outcm1("outcm1.dat");
	ofstream outcm2("outcm2.dat");
   //montagem de conectividades entre elementos
   firstmesh->BuildConnectivity();
   //malha computacional
   TPZCompMesh *secondmesh = new TPZCompMesh(firstmesh);
   secondmesh->SetName("Malha Computacional : Conectividades e Elementos");
   //material
   TPZMaterial *pl = LerMaterial("placa1.dat");
   secondmesh->InsertMaterialObject(pl);
   //forcingfunction
   //CC : condicões de contorno
   TPZBndCond *bc;
   REAL big = 1.e12;
   TPZFMatrix val1(6,6,0.),val2(6,1,0.);
   val1(0,0)=big;
   val1(1,1)=big;
   val1(2,2)=big;
   val1(3,3)=big;
   val1(4,4)=big;
   val1(5,5)=big;

 	TPZGeoElBC(elg0,4,-1,*firstmesh);
   bc = pl->CreateBC(-1,0,val1,val2);
   secondmesh->InsertMaterialObject(bc);

   val1.Zero();
	TPZGeoElBC(elg0,5,-2,*firstmesh);
   bc = pl->CreateBC(-2,2,val1,val2);
   secondmesh->InsertMaterialObject(bc);

 	TPZGeoElBC(elg0,6,-3,*firstmesh);
   bc = pl->CreateBC(-3,2,val1,val2);
   secondmesh->InsertMaterialObject(bc);

   val1.Zero();
	TPZGeoElBC(elg0,7,-4,*firstmesh);
   bc = pl->CreateBC(-4,2,val1,val2);
   secondmesh->InsertMaterialObject(bc);

   //ordem de interpolacao
   int ord;
   cout << "Entre ordem 1,2,3,4,5 : ";
   cin >> ord;
//   TPZCompEl::gOrder = ord;
   cmesh.SetDefaultOrder(ord);
   //construção malha computacional
   secondmesh->AutoBuild();
   //redistribuicao de ordem aos lados do elemento
	int nel = secondmesh->ElementVec().NElements();
   cout << "\nEntre Ordem por Elemento\n";
   for(int index = 0;index<nel;index++) {
      TPZInterpolatedElement *el;
		el = ((TPZInterpolatedElement *) secondmesh->ElementVec()[index]);
      cout << "\nElemento Geometrico de Id " << el->Reference()->Id() << " : ordem -> ";
      cin >> ord;
      if(el) el->PRefine(ord);
   }
   TPZVec<int> csub(0);
   int n1=1,n2;
   while(n1) {
	   cout << "Id do elemento geometrico a dividir ? : ";
      cin >> n1;
      if(n1 < 0) break;
      int nelc = secondmesh->ElementVec().NElements();
      int el;
      TPZCompEl *cpel;
      for(el=0;el<nelc;el++) {
         cpel = secondmesh->ElementVec()[el];
         if(cpel && cpel->Reference()->Id() == n1) break;
      }
      secondmesh->Divide(el,csub,1);
      n1 = 1;
   }
/*   cout << "Entre numeros minimos e maximos de elementos (20-50) : ";
   cin >>n1 >>n2;
   //numeros minimos e maximos de elementos que a malha nao devera ultrapassar apos agrupamento e divisao, respectivamente
   CycleRefinements(*secondmesh,n1,n2);*/
   //analysis
	secondmesh->InitializeBlock();
   secondmesh->Print(outcm1);
   TPZAnalysis an(secondmesh,outcm1);
   int numeq = secondmesh->NEquations();
   secondmesh->Print(outcm1);
   outcm1.flush();
/*   TPZVec<int> skyline;
   secondmesh->Skyline(skyline);
	TPZSkylMatrix *stiff = new TPZSkylMatrix(numeq,skyline);*/
   TPZFMatrix *stiff = new TPZFMatrix(numeq,numeq);
   an.SetMatrix(stiff);
   an.Solver().SetDirect(ELU);
   //an.Solver().SetDirect(ELDLt);
   //an.Solver().SetDirect(ECholesky);
   //an.Solver().SetJacobi(4, 1E-8, 0);
   //an.Solver().SetSOR(4, 0.2, 1E-8, 0);
   //an.Solver().SetSSOR(6, 1.3, 1E-8,0);
   secondmesh->SetName("Malha Computacional :  Connects e Elementos");
   // Posprocessamento
   an.Run(outcm2);
   TPZVec<char *> scalnames(3);
   scalnames[0] = "Deslocx";
   scalnames[1] = "Deslocy";
   scalnames[2] = "Deslocz";
   TPZVec<char *> vecnames(0);
   char plotfile[] =  "plot.pos";
   an.DefineGraphMesh(2, scalnames, vecnames, plotfile);
   an.Print("FEM SOLUTION ",outcm1);
   firstmesh->Print(outgm1);
   outgm1.flush();
   delete secondmesh;
   delete firstmesh;
   return 0;
}
//FIM MAIN FIM MAIN FIM MAIN FIM MAIN FIM MAIN FIM MAIn FIM MAIN FIM MAIN
TPZMaterial *LerMaterial(char *filename) {
	ifstream input(filename);
   TPZFMatrix naxes(3,3);
   REAL ni1,ni2,h,E1,E2,G12,G13,G23,f;
   REAL n00,n01,n02,n10,n11,n12,n20,n21,n22;
   TPZVec<REAL> xf(6);
   int matindex;
   input >> matindex;
   input >> f  >>  h  >>
   		  E1  >>  E2 >>
           G12 >> G13 >> G23 >>
           ni1 >> ni2;
   input >> n00 >> n01 >> n02 >> n10 >> n11 >> n12 >> n20 >> n21 >> n22;
	input >> xf[0] >> xf[1] >> xf[2] >> xf[3] >> xf[4] >> xf[5];
   naxes(0,0) =  n00;    naxes(0,1) =  n01;    naxes(0,2) =  n02;
   naxes(1,0) =  n10;    naxes(1,1) =  n11;    naxes(1,2) =  n12;
   naxes(2,0) =  n20;    naxes(2,1) =  n21;    naxes(2,2) =  n22;
   return new TPZPlaca(matindex,h,f,E1,E2,ni1,ni2,G12,G13,G23,naxes,xf);
}


/*
TESTE 1 TESTE 1 TESTE 1 TESTE 1 TESTE 1 TESTE 1 TESTE 1

   val1(0,0)=big;
   val1(1,1)=big;
   val1(2,2)=big;
   val1(3,3)=big;
   val1(4,4)=big;
   val1(5,5)=big;
 	TPZGeoElBC(elg0,7,-2,*firstmesh);
   bc = pl->CreateBC(-2,2,val1,val2);
   secondmesh->InsertMaterialObject(bc);

// 	TPZGeoElBC(elg0,0,-3,*firstmesh);//equivale a val1(1,1)=big;
//   bc = pl->CreateBC(-3,0,val1,val2);
//   secondmesh->InsertMaterialObject(bc);

   val1.Zero();
   val2(0,0) = 100.;
	TPZGeoElBC(elg0,5,-1,*firstmesh);
   bc = pl->CreateBC(-1,1,val1,val2);
   secondmesh->InsertMaterialObject(bc);
*/

void ExactSol(TPZVec<REAL> &loc,TPZVec<REAL> &val,TPZMatrix &deriv) {

   val[0]  = 4.*loc[0];//E = 100.
   val[1]  = -loc[1];//-.25*loc[1];
   deriv(0,0) = 4.;
   deriv(1,0) = 0.;
   deriv(0,1) = 0.;
   deriv(1,1) = -1.;//-.25
}

