

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

#include "pzelgq2d.h"//Cedric

#define NOTDEBUG
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


int main() {
   //malha geometrica
   TPZGeoMesh *firstmesh = new TPZGeoMesh;
   firstmesh->NodeVec().Resize(4);
   TPZVec<REAL> coord(2);
   coord[0] = 0.;
   coord[1] = 0.;
   //nos geometricos
   firstmesh->NodeVec()[0].Initialize(coord,*firstmesh);
   coord[0] = 1.0;
   firstmesh->NodeVec()[1].Initialize(coord,*firstmesh);
   coord[1] = 1.0;
   firstmesh->NodeVec()[2].Initialize(coord,*firstmesh);
   coord[0] = 0.0;
   firstmesh->NodeVec()[3].Initialize(coord,*firstmesh);
   TPZVec<int> nodeindexes(4);
   nodeindexes[0] = 0;//local[i] = global[i] , i=0,1,2,3
   nodeindexes[1] = 1;
   nodeindexes[2] = 2;
   nodeindexes[3] = 3;
   //elementos geometricos
   TPZGeoElQ2d *elq1 = new TPZGeoElQ2d(nodeindexes,1,*firstmesh);
/*   int i,sen;;
   cout<<"Entre i = 0,1,2,3 ";
   cin>>i;
   cout<<"Sentido local direito/inverso : 0/1 ?  ";
   cin>>sen;
   if(sen==0) {//direito
        nodeindexes[0] = (0+i)%4;//local[i] = global[j] , i,j em {0,1,2,3}
        nodeindexes[1] = (1+i)%4;
        nodeindexes[2] = (2+i)%4;
        nodeindexes[3] = (3+i)%4;
	} else {//inverso
        nodeindexes[0] = (3+i)%4;//local[i] = global[j] , i,j em {0,1,2,3}
        nodeindexes[1] = (2+i)%4;
        nodeindexes[2] = (1+i)%4;
        nodeindexes[3] = (0+i)%4;
   }
   TPZGeoElQ2d *elq2 = new TPZGeoElQ2d(nodeindexes,1,*firstmesh);   */
   //Arquivos de saida
	ofstream outgm1("outgm1.dat");
   ofstream outcm1("outcm1.dat");
   //montagem de conectividades entre elementos
   firstmesh->BuildConnectivity();
  	//teste de divisao geometrica
/*   TPZVec<TPZGeoEl *> vecsub,vecsub1;
   elq1->Divide(vecsub);//divide 0
   vecsub[0]->Divide(vecsub1);//divide 2
	elq2->Divide(vecsub1);//div 1
   cout<<"Entre filho a ser subdividido  ";
   cin>>i;
   vecsub1[i]->Divide(vecsub1);//div filho i do 2*/
 	firstmesh->Print(outgm1);
   outgm1.flush();
   delete firstmesh;
   return 0;
}
//FIM MAIN FIM MAIN FIM MAIN FIM MAIN FIM MAIN FIM MAIN FIM MAIN FIM MAIN
