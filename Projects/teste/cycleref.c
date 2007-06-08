

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
int Coarsen(TPZCompEl *cel,TPZCompMesh &mesh);
void DivideAny(TPZCompMesh &cmsh);
void CoarsenAny(TPZCompMesh &cmsh);
void CycleRefinements(TPZCompMesh &cmesh, int minel, int maxel);
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
   coord[0] = 1.0;//10
   firstmesh->NodeVec()[1].Initialize(coord,*firstmesh);
   coord[0] = 1.0;
   coord[1] = 1.0;//11
   firstmesh->NodeVec()[2].Initialize(coord,*firstmesh);
   coord[0] = 0.0;//01
   firstmesh->NodeVec()[3].Initialize(coord,*firstmesh);
   TPZVec<int> nodeindexes(3);//triangulo
   nodeindexes[0] = 0;//local[i] = global[i] , i=0,1,2,3
   nodeindexes[1] = 1;
   nodeindexes[2] = 2;//ex. anterior ; 021
   //elementos geometricos
   TPZGeoElT2d *elg0 = new TPZGeoElT2d(nodeindexes,1,*firstmesh);
   //orientacao local de um segundo elemento superposto
   nodeindexes[0] = 0;
   nodeindexes[1] = 2;
   nodeindexes[2] = 3;//302
   TPZGeoElT2d *elg1 = new TPZGeoElT2d(nodeindexes,1,*firstmesh);
   //Arquivos de saida
	ofstream outgm1("outgm1.dat");
   ofstream outcm1("outcm1.dat");
	ofstream outcm2("outcm2.dat");
   //montagem de conectividades entre elementos
   firstmesh->BuildConnectivity();
  	//teste de divisao geometrica
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
   //CC : condicao de contorno
 /*	TPZGeoElBC(elg1,4,-1,*firstmesh);
   TPZFMatrix val1(1,1,0.), val2(1,1,1.);
   int bcindex;
   TPZBndCond *bc = mat->CreateBC(-1,0,val1,val2);
   bcindex = secondmesh->InsertMaterialObject(bc);
   if(carga==0) bc->SetForcingFunction(force);
   else bc->SetForcingFunction(bcforce);  */
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
   int intinpose[] = {2,2,2,2,2,2,2};//{4,4,4,4,4,4,4};
//   int randseed;
//   cout << "enter rand seed\n";
//   cin >> randseed;
//   srand(randseed);
	randomize();
   for(int index=0;index<nel;index++) {
   	//int ran = random(3)+2;
      TPZInterpolatedElement *el=0;
      el = ((TPZInterpolatedElement *) secondmesh->ElementVec()[index]);
      el->PRefine(intinpose[index]);
      //el->PRefine(ran);
   }
   TPZVec<int> csub;
   int n1=1,n2;
   while(n1) {
	   cout << "Id geometrico do elemento a dividir ? : ";
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
   TPZVec<char *> scalnames(1);
   scalnames[0] = "state";
   TPZVec<char *> vecnames(0);
   char plotfile[] =  "plot.plt";
   an.DefineGraphMesh(2, scalnames, vecnames, plotfile);
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
int Coarsen(TPZCompEl *cel,TPZCompMesh &mesh) {
	if(!cel) return 0;
	TPZGeoEl *ref = cel->Reference();
   TPZGeoEl *father = ref->Father();
   if(!father) return 0;
   int numsub = father->NSubElements();
   TPZVec<TPZCompEl *> compsubel(numsub);
   int returncode = 1;
   for(int is=0; is<numsub; is++) {
		TPZGeoEl *gsub = father->SubElement(is);
      compsubel[is] = gsub->Reference();
      if(!gsub->Reference()) returncode = 0;

   }
   if(!returncode) return returncode;
   TPZManVector<int> subindex;
   int index;
   numsub = compsubel.NElements();
   subindex.Resize(numsub);
   for(int is=0; is<numsub; is++) subindex[is] = compsubel[is]->Index();
   if(returncode) {
   	cout << "Coarsening " << compsubel[0]->Reference()->Id() << ' ' << compsubel[1]->Reference()->Id() << endl;
   	mesh.Coarsen(subindex,index);
      cout << "Created " << mesh.ElementVec()[index]->Reference()->Id() << endl;
   }
   return returncode;
}

ofstream locfile("CycleRefinements");

void DivideAny(TPZCompMesh &cmesh) {
   TPZAdmChunkVector<TPZCompEl *> &vec = cmesh.ElementVec();
   int nelem = vec.NElements();
	int elindex = rand()%nelem;
   TPZCompEl *cel = vec[elindex];
   while(!cel) {
   	elindex = rand()%nelem;
   	cel = vec[elindex];
   }
   TPZManVector<int> subindex;
   cout << "Dividing " << vec[elindex]->Reference()->Id() << endl;
   cout.flush();
   if(vec[elindex]->Reference()->Id() == 4) {
	   cmesh.Reference()->Print(locfile);
	   cmesh.Print(locfile);
   }

   cmesh.Divide(elindex,subindex);
   cout << " creating " << vec[subindex[0]]->Reference()->Id() << ' ' << vec[subindex[1]]->Reference()->Id() << endl;
   //cmesh.Reference()->Print(locfile);
//   locfile.flush();
//   cmesh.CleanUpUnconnectedNodes();
}

void CoarsenAny(TPZCompMesh &cmesh) {
   TPZAdmChunkVector<TPZCompEl *> &vec = cmesh.ElementVec();
   int nelem = vec.NElements();
	int elindex = rand()%nelem;
   TPZCompEl *cel = vec[elindex];
   while(!cel || !Coarsen(cel,cmesh)) {
		elindex = rand()%nelem;
   	cel = vec[elindex];
   }
//   cmesh.ExpandSolution();
//   cmesh.CleanUpUnconnectedNodes();
}

ofstream logfile("log.dat");

void CycleRefinements(TPZCompMesh &cmesh, int minel, int maxel) {
   int cycle;
   int nelem = cmesh.NElements();
   int i,totalel;
   TPZAdmChunkVector<TPZCompEl *> &vec = cmesh.ElementVec();
	for(cycle = 0; cycle<10; cycle++) {
   	while(nelem < maxel) {
      	DivideAny(cmesh);
         nelem=0;
         totalel = vec.NElements();
         for(i=0; i<totalel; i++) if(vec[i]) nelem++;
      }

      while(nelem > minel) {
      	CoarsenAny(cmesh);
         nelem=0;
         totalel = vec.NElements();
         for(i=0; i<totalel; i++) if(vec[i]) nelem++;
      }
      cmesh.ComputeNodElCon();
//      cmesh.Print(logfile);
//      logfile.flush();
      cmesh.ExpandSolution();
      cmesh.CleanUpUnconnectedNodes();
      cmesh.ExpandSolution();
//      cmesh.Print(logfile);
//      logfile.flush();
//      cmesh.ComputeConnectSequence();
      if(cycle == 0) cmesh.Print(logfile);
      logfile.flush();
   }
//   cmesh.Print(logfile);
//   logfile.flush();
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
