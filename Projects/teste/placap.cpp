

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

#include "pzelgq2d.h"
#include "pzelgt2d.h"
#include "pzmat2dlin.h"
#include "pzanalysis.h"
#include "pzmetis.h"

#include <stdio.h>
#include <time.h>
#include "pzelct2d.h"

#define NOTDEBUG
int Coarsen(TPZCompEl *cel,TPZCompMesh &mesh);
void DivideAny(TPZCompMesh &cmsh);
void CoarsenAny(TPZCompMesh &cmsh);
void CycleRefinements(TPZCompMesh &cmesh, int minel, int maxel);
void force(TPZVec<REAL> &x, TPZVec<REAL> &f);
void bcforce(TPZVec<REAL> &x, TPZVec<REAL> &f);
void derivforce(TPZVec<REAL> &x, TPZVec<REAL> &f);
void ExactSol(TPZVec<REAL> &loc,TPZVec<REAL> &val,TPZMatrix &deriv);
int main() {

   //malha geometrica
   TPZGeoMesh *firstmesh = new TPZGeoMesh;
   firstmesh->SetName("Malha Geometrica : Nós e Elementos");
   firstmesh->NodeVec().Resize(3);
   TPZVec<REAL> coord(2);
   coord[0] = 0.;
   coord[1] = 0.;
   //nos geometricos
   firstmesh->NodeVec()[0].Initialize(coord,*firstmesh);
   coord[0] = 1.0;
   firstmesh->NodeVec()[1].Initialize(coord,*firstmesh);
   coord[0] = 0.0;
   coord[1] = 1.0;
   firstmesh->NodeVec()[2].Initialize(coord,*firstmesh);
   coord[0] = 0.0;
   //firstmesh->NodeVec()[3].Initialize(coord,*firstmesh);
   //TPZVec<int> nodeindexes(3);
   nodeindexes[0] = 1;
   nodeindexes[1] = 2;
   nodeindexes[2] = 0;
   //elementos geometricos
   TPZGeoElT2d *elg0 = new TPZGeoElT2d(nodeindexes,1,*firstmesh);
///////////////////////////////////////////////////////
   //nodeindexes[0] = 1;
   //nodeindexes[1] = 2;
   //nodeindexes[2] = 0;
   //TPZGeoElT2d *elg0 = new TPZGeoElT2d(nodeindexes,1,*firstmesh);
///////////////////////////////////////////////////////
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
   TPZMaterial *mat = new TPZElasticityMaterial(1, 25.,.25,0.,0.);
   secondmesh->InsertMaterialObject(mat);
   //forcingfunction
   //CC : condicões de contorno
   TPZBndCond *bc;
   REAL big = 1.e12;
   TPZFMatrix val1(2,2,0.), val2(2,1,0.);

   val1(1,1) = big;
 	TPZGeoElBC(elg0,5,-1,*firstmesh);
   bc = mat->CreateBC(-1,2,val1,val2);
   secondmesh->InsertMaterialObject(bc);

   val1(1,1) = 0.;
   val2(0,0) =  70.71067811865;
	TPZGeoElBC(elg0,3,-2,*firstmesh);
   bc = mat->CreateBC(-2,1,val1,val2);
   secondmesh->InsertMaterialObject(bc);

   val1(0,0) = big;
   val2(0,0) = .0;
   TPZGeoElBC(elg0,4,-3,*firstmesh);
   bc = mat->CreateBC(-3,2,val1,val2);
   secondmesh->InsertMaterialObject(bc);

   //ordem de interpolacao
   int ord;
   cout << "Entre ordem 1,2,3,4,5 : ";
   cin >> ord;
   //TPZCompEl::gOrder = ord;
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
   TPZVec<int> csub;
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
   TPZVec<char *> scalnames(1);
   scalnames[0] = "state";
   TPZVec<char *> vecnames(0);
   char plotfile[] =  "plot.plt";
   an.DefineGraphMesh(2, scalnames, vecnames, plotfile);
   an.Print("FEM SOLUTION ",outcm1);
   firstmesh->Print(outgm1);
   outgm1.flush();
///////////////////////////////////////////////////////
   REAL estimate,true_error,L2_error;
   TPZVec<REAL> flux(0);
   secondmesh->EvaluateError(ExactSol,true_error,L2_error,estimate);
   ofstream error("Error.dat");
   error << "True Error    : " << true_error << endl;
   error << "L2 Error      : " << L2_error << endl;
   error << "estimate Error : " << estimate << endl;
   error.flush();
///////////////////////////////////////////////////////
   delete secondmesh;
   delete firstmesh;
   return 0;
}
//FIM MAIN FIM MAIN FIM MAIN FIM MAIN FIM MAIN FIM MAIn FIM MAIN FIM MAIN
void ExactSol(TPZVec<REAL> &loc,TPZVec<REAL> &val,TPZMatrix &deriv) {

   val[0]  = 4.*loc[0];//E = 100.
   val[1]  = -loc[1];//-.25*loc[1];
   deriv(0,0) = 4.;
   deriv(1,0) = 0.;
   deriv(0,1) = 0.;
   deriv(1,1) = -1.;//-.25
}

void force(TPZVec<REAL> &x, TPZVec<REAL> &f) {
   f[0]  = 0.;
}

void bcforce(TPZVec<REAL> &x, TPZVec<REAL> &f) {
	f[0] = 0.;
}
void derivforce(TPZVec<REAL> &x, TPZVec<REAL> &f) {

	cout << "\nmain()::derivforce : Not implemented yet\n";
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

/*
   val1(1,1) = big;
 	TPZGeoElBC(elg0,5,-1,*firstmesh);
   bc = mat->CreateBC(-1,2,val1,val2);
   secondmesh->InsertMaterialObject(bc);

   val1(1,1) = 0.;
   val2(0,0) = 100.0;
	TPZGeoElBC(elg0,3,-2,*firstmesh);
   bc = mat->CreateBC(-2,1,val1,val2);
   secondmesh->InsertMaterialObject(bc);

   val1(0,0) = big;
   val2(0,0) = .0;
   TPZGeoElBC(elg1,3,-3,*firstmesh);
   bc = mat->CreateBC(-3,2,val1,val2);
   secondmesh->InsertMaterialObject(bc);
*/
