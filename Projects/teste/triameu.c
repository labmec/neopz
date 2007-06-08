

#include "pzfmatrix.h"
#include "pzskylmat.h"
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
#include "pzanalysis.h"
#include "pzmetis.h"

#include <stdio.h>//Cedric
#include <time.h>//Cedric
#include "pzelct2d.h"//Cedric

#define NOTDEBUG
void force(TPZVec<REAL> &x, TPZVec<REAL> &f);
void bcforce(TPZVec<REAL> &x, TPZVec<REAL> &f);
void derivforce(TPZVec<REAL> &x, TPZVec<REAL> &f);

int main() {

   //malha geometrica
   TPZGeoMesh *firstmesh = new TPZGeoMesh;
   firstmesh->SetName("Malha Geometrica : Nós e Elementos");
   firstmesh->NodeVec().Resize(4);
   TPZVec<REAL> coord(2);
   coord[0] = 0.;
   coord[1] = 0.;//00
   //nos geometricos
   firstmesh->NodeVec()[0].Initialize(coord,*firstmesh);
   coord[0] = 2.0;//10
   firstmesh->NodeVec()[1].Initialize(coord,*firstmesh);
   coord[0] = 2.0;
   coord[1] = 2.0;//11
   firstmesh->NodeVec()[2].Initialize(coord,*firstmesh);
   coord[0] = 0.0;//01
   firstmesh->NodeVec()[3].Initialize(coord,*firstmesh);
   TPZVec<int> nodeindexes(3);//triangulo
   nodeindexes[0] = 0;//local[i] = global[i] , i=0,1,2,3
   nodeindexes[1] = 1;
   nodeindexes[2] = 3;//ex. anterior ; 021
   //elementos geometricos
   TPZGeoElT2d *elg0 = new TPZGeoElT2d(nodeindexes,1,*firstmesh);
   //orientacao local de um segundo elemento superposto
   nodeindexes[0] = 2;
   nodeindexes[1] = 3;
   nodeindexes[2] = 1;//302
   TPZGeoElT2d *elg1 = new TPZGeoElT2d(nodeindexes,1,*firstmesh);
   nodeindexes[0] = 1;
   nodeindexes[1] = 0;
   nodeindexes[2] = 2;//102
//   TPZGeoElT2d *elg2 = new TPZGeoElT2d(nodeindexes,1,*firstmesh);
   nodeindexes[0] = 2;
   nodeindexes[1] = 0;
   nodeindexes[2] = 3;//203
//   TPZGeoElT2d *elg3 = new TPZGeoElT2d(nodeindexes,1,*firstmesh);
   //Arquivos de saida
	ofstream outgm1("outgm1.dat");
   ofstream outcm1("outcm1.dat");
	ofstream outcm2("outcm2.dat");
   //montagem de conectividades entre elementos
   firstmesh->BuildConnectivity();
  	//teste de divisao geometrica : 1 elemento
   TPZVec<TPZGeoEl *> vecsub,vecsub1;
   elg0->Divide(vecsub);
   //1o exemplo
//   elg2->Divide(vecsub);//2:4567
//   elg1->Divide(vecsub1);//1:891011
//   vecsub[2]->Divide(vecsub);//6:12131415
//	vecsub[1]->Divide(vecsub);//13:16171819
//   elg0->Divide(vecsub);//0:20212223
//   vecsub[3]->Divide(vecsub);//23:24252627
	//2o exemplo
//   elg2->Divide(vecsub);//2:4567
//   elg1->Divide(vecsub1);//1:891011
//   vecsub[0]->Divide(vecsub);//4:12131415
//	vecsub[1]->Divide(vecsub);//13:16171819
//   elg0->Divide(vecsub);//0:20212223
//   vecsub[3]->Divide(vecsub);//23:24252627
   //malha computacional
   TPZCompMesh *secondmesh = new TPZCompMesh(firstmesh);
   secondmesh->SetName("Malha Computacional : Conectividades e Elementos");
   //material
   int matindex;
   TPZFMatrix k(1,1,1.),f(1,1,0.),c(1,2,0.);
   TPZMat2dLin * mat = new TPZMat2dLin(1);
   matindex = secondmesh->InsertMaterialObject(mat);
   mat->SetMaterial(k,c,f);
   //forcingfunction
   int carga;
   cout << "Entre teste funcao/derivada (0/1) : ";
   cin >> carga;
   if(carga==0) {
   	mat->SetForcingFunction(force);
   } else {//default
   	mat->Ck()(0,0) = 1;
   	mat->Ck()(0,1) = 1;
   	mat->SetForcingFunction(derivforce);

   }
   cout << "CC ? (0/1) : ";
   int resp;
   cin >> resp;
   if(resp==1) {
        //CC : condicao de contorno
        TPZGeoElBC(vecsub[0],5,-1,*firstmesh);
        TPZGeoElBC(vecsub[2],5,-1,*firstmesh);
        TPZFMatrix val1(1,1,0.), val2(1,1,1.);
        int bcindex;
        TPZBndCond *bc = mat->CreateBC(-1,0,val1,val2);
        bcindex = secondmesh->InsertMaterialObject(bc);
        if(carga==0) bc->SetForcingFunction(bcforce);
        else bc->SetForcingFunction(bcforce);
   }
   //ordem de interpolacao
   int ord;
   cout << "Entre ordem 1,2,3,4,5 : ";
   cin >> ord;
//   TPZCompEl::gOrder = ord;
   cmesh.SetDefaultOrder(ord);
   //construção malha computacional
   secondmesh->AutoBuild();
   //redistribuicao de ordem aos lados do elemento
	int index,nel = secondmesh->ElementVec().NElements();
   cout << "Ordem aleatoria/imposta (1/0) ? ";
   ord;
   cin >> ord;
   if(ord==1) {
        randomize();
        for(index=0;index<nel;index++) {
           int ran = random(3)+2;
           TPZInterpolatedElement *el=0;
           el = ((TPZInterpolatedElement *) secondmesh->ElementVec()[index]);
           el->PRefine(ran);
        }
	} else {
      //int *intinpose = new int[nel];
      cout << "\nEntre Ordem por Elemento : 1,2,3,4,5\n";
      for(index = 0;index<nel;index++) {
      	cout << "\nIndex do Elemento " << index << " ordem -> ";
		   cin >> ord;//intinpose[index];
         TPZInterpolatedElement *el=0;
         el = ((TPZInterpolatedElement *) secondmesh->ElementVec()[index]);
         el->PRefine(ord);
      }
   }
   secondmesh->Print(outcm1);
   outcm1.flush();
   //analysis
	secondmesh->InitializeBlock();
   TPZAnalysis an(secondmesh);
   int numeq = secondmesh->NEquations();
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
   an.Print("FEM SOLUTION ",outcm1);
   firstmesh->Print(outgm1);
   outgm1.flush();
///////////////////////////////////////////////////////
   delete secondmesh;
   delete firstmesh;
   return 0;
}




//FIM MAIN FIM MAIN FIM MAIN FIM MAIN FIM MAIN FIM MAIN FIM MAIN FIM MAIN
void bcforce(TPZVec<REAL> &x, TPZVec<REAL> &f) {
	f[0] = 0.0;
}
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



/*   TPZStack<int> elgraph(0);
   TPZStack<int> elgraphindex(0);
   int nnod;
   secondmesh->ComputeElGraph(elgraph,elgraphindex,nnod);
   int nel = elgraphindex.NElements()-1;
   TPZMetis metis(nel,nnod);
   metis.SetElementGraph(elgraph,elgraphindex);
   TPZVec<int> perm(0), iperm(0);
   metis.Resequence(perm,iperm);
   secondmesh->Permute(perm);
*/

/*   TPZVec<char *> scalnames(3),vecnames(1);
   scalnames[0] = "SigmaX";
   scalnames[1] = "SigmaY";
   scalnames[2] = "TauXY";
   vecnames[0]  = "Displacement";
   //vecnames[1]  = "PrincipalStresses";
   an.DefineGraphMesh(2,scalnames,vecnames,"MvwModel.pos");
*/

/*   int bcindex = firstmesh->BCElementVec().AllocateNewElement();
   firstmesh->BCElementVec()[bcindex] = bcel1;//elemento geometrico
   TPZBndCond *bc = mat->CreateBC(-1,0,val1,val2);
   bcindex = secondmesh->BndCondVec().AllocateNewElement();
   secondmesh->BndCondVec()[bcindex] = bc;//BC : AutoBuidl gerara o elemento computacional

   bcindex = firstmesh->BCElementVec().AllocateNewElement();
   firstmesh->BCElementVec()[bcindex] = bcel2;
   bc = mat->CreateBC(-1,0,val1,val2);
   bcindex = secondmesh->BndCondVec().AllocateNewElement();
   secondmesh->BndCondVec()[bcindex] = bc;
*/


