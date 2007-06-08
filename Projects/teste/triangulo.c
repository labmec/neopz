

#include "pzfmatrix.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzgnode.h"
#include "pzsolve.h"
#include "pzelg1d.h"
#include "pzmat1dlin.h"
#include "pzelmat.h"
#include "pzcompel.h"
#include "pzbndcond.h"
#include "pzstack.h"
#include "pztrnsform.h"
#include <stdlib.h>
#include <iostream.h>

#include "pzelgq2d.h"
#include "pzelgt2d.h"
#include "pzmat2dlin.h"

#define NOTDEBUG
void force(TPZVec<REAL> &x, TPZVec<REAL> &f);
void derivforce(TPZVec<REAL> &x, TPZVec<REAL> &f);

int main() {

   //malha geometrica
   TPZGeoMesh *firstmesh = new TPZGeoMesh;
   firstmesh->NodeVec().Resize(4);
   TPZVec<REAL> coord(2);
   coord[0] = 0.;
   coord[1] = 0.;//00
   //nos geometricos
   firstmesh->NodeVec()[0].Initialize(coord,*firstmesh);
   coord[0] = 1.0;//10
   firstmesh->NodeVec()[1].Initialize(coord,*firstmesh);
   coord[0] = 1.0;
   coord[1] = 1.0;//11
   firstmesh->NodeVec()[2].Initialize(coord,*firstmesh);
   coord[0] = 0.0;//01
   firstmesh->NodeVec()[3].Initialize(coord,*firstmesh);
   TPZVec<int> nodeindexes(3);//triangulo
   nodeindexes[0] = 0;//local[i] = global[i] , i=0,1,2,3
   nodeindexes[1] = 2;
   nodeindexes[2] = 1;
   //elementos geometricos
   TPZGeoElT2d *elg0 = new TPZGeoElT2d(nodeindexes,1,*firstmesh);
   //orientacao local de um segundo elemento superposto
   nodeindexes[0] = 3;
   nodeindexes[1] = 0;
   nodeindexes[2] = 2;
   TPZGeoElT2d *elg1 = new TPZGeoElT2d(nodeindexes,1,*firstmesh);
   nodeindexes[0] = 1;
   nodeindexes[1] = 0;
   nodeindexes[2] = 2;
   TPZGeoElT2d *elg2 = new TPZGeoElT2d(nodeindexes,1,*firstmesh);
   nodeindexes[0] = 2;
   nodeindexes[1] = 0;
   nodeindexes[2] = 3;
   TPZGeoElT2d *elg3 = new TPZGeoElT2d(nodeindexes,1,*firstmesh);
   //Arquivos de saida
	ofstream outgm1("outgm1.dat");
   ofstream outcm1("outcm1.dat");
	ofstream outcm2("outcm2.dat");
   //montagem de conectividades entre elementos
   firstmesh->BuildConnectivity();
 	firstmesh->Print(outgm1);
   outgm1.flush();
  	//teste de divisao geometrica : 1 elemento
   TPZVec<TPZGeoEl *> vecsub,vecsub1;
 	firstmesh->Print(outgm1);
   outgm1.flush();
   elg2->Divide(vecsub);//2:4567
   elg1->Divide(vecsub1);//1:891011
   vecsub[2]->Divide(vecsub);//6:12131415
	vecsub[1]->Divide(vecsub);//13:16171819
   elg0->Divide(vecsub);//0:20212223
   vecsub[3]->Divide(vecsub);//23:24252627
//   vecsub[0]->Divide(vecsub);//24:28293031
 	firstmesh->Print(outgm1);
   outgm1.flush();
   //malha computacional
   TPZCompMesh *secondmesh = new TPZCompMesh(firstmesh);
   //material
   int matindex = secondmesh->MaterialVec().AllocateNewElement();
   TPZFMatrix k(1,1,1.),f(1,1,0.),c(1,2,0.);
   TPZMat2dLin * mat = new TPZMat2dLin(1);
   mat->SetMaterial(k,c,f);
   //forcingfunction
   int carga;
   cout << "Entre teste funcao/derivada (0/1) : ";
   cin >> carga;
   if(carga==0) mat->SetForcingFunction(force);
   if(carga==1) {
   	mat->Ck()(0,0) = 1;
   	mat->Ck()(0,1) = 1;
   	mat->SetForcingFunction(derivforce);

   }
   secondmesh->MaterialVec()[matindex] = mat;
   //CC : condicao de contorno
   //ordem de interpolacao
   int ord;
   cout << "Entre ordem 1,2,3,4,5 : ";
   cin >> ord;
//   TPZCompEl::gOrder = ord;
   cmesh.SetDefaultOrder(ord);
   //constroe a malha computacional
   secondmesh->AutoBuild();
   secondmesh->Print(outcm1);
   outcm1.flush();
   secondmesh->InitializeBlock();
   secondmesh->ComputeConnectSequence();
	//Resolucao do sistema
   TPZFMatrix Rhs(secondmesh->NEquations(),1),Stiff(secondmesh->NEquations(),secondmesh->NEquations()),U;
   Stiff.Zero();
   Rhs.Zero();
   secondmesh->Assemble(Stiff,Rhs);
   Rhs.Print("Rhs teste",outcm2);
   Stiff.Print("Bloco teste",outcm2);
	Rhs.Print("Computational Mesh -> fBlock",outcm2);
   TPZMatrixSolver solver(&Stiff);
   solver.SetDirect(ELU);
   solver.Solve(Rhs,U);
   U.Print("Resultado",outcm2);
   secondmesh->LoadSolution(U);
   secondmesh->Solution().Print("Mesh solution ",outcm2);
   secondmesh->ConnectSolution(outcm2);
   //   TPZElementMatrix ek,ef;
//   secondmesh->ElementVec()[0]->CalcStiff(ek,ef);
//	ek.fMat->Print();
//   ef.fMat->Print();
   delete secondmesh;
   delete firstmesh;
   return 0;
}
//FIM MAIN FIM MAIN FIM MAIN FIM MAIN FIM MAIN FIM MAIN FIM MAIN FIM MAIN
void force(TPZVec<REAL> &x, TPZVec<REAL> &f) {
	static int init=0,n,m;
   if(init == 0) {
  		init = 1;
      cout<<"\nMaterial::TPZMat2dLin -> Teste com monomios xn*ym. Entre expoentes n,m\n";
      cin>>n>>m;
	}
   int r = f.NElements();
   int ic;
   for(ic=0; ic< r; ic++) {
      if (!n && !m) f[ic]  = 1.;
      if ( n      ) f[ic]  = pow(x[0], n);
      if ( n &&  m) f[ic] *= pow(x[1], m);
      if (!n &&  m) f[ic]  = pow(x[1], m);
   }
}

void derivforce(TPZVec<REAL> &x, TPZVec<REAL> &f) {
	static int init=0,n,m;
   if(init == 0) {
  		init = 1;
      cout<<"\nMaterial::TPZMat2dLin -> Solucao  f =  xn*ym. Entre expoentes nao negativos n,m \n";
      cin>>n>>m;
	}
   int r = f.NElements();
   int ic;
   for(ic=0; ic< r; ic++) {
   	if(n>0 && m>0) {
      	f[ic]  = pow(x[0], n) * pow(x[1], m);
      	f[ic] += n*pow(x[0], n-1)*pow(x[1], m) + m*pow(x[0], n)*pow(x[1], m-1);
      } else
   	if(n==0 && m>0) {
      	f[ic]  = pow(x[1], m);
      	f[ic] += m*pow(x[1], m-1);
      } else
   	if(n>0 && m==0) {
      	f[ic]  = pow(x[0], n);
      	f[ic] += n*pow(x[0], n-1);
      } else {//outros casos
      	f[ic] = 0;
      }
   }
}








