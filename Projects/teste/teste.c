
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

#define NOTDEBUG

int Coarsen(TPZCompEl *cel,TPZCompMesh &mesh);
void DivideAny(TPZCompMesh &cmsh);
void CoarsenAny(TPZCompMesh &cmsh);
void CycleRefinements(TPZCompMesh &cmesh, int minel, int maxel);

int main() { // int main(int argc,char *argv[]) {
   //malha geometrica
   TPZGeoMesh *firstmesh = new TPZGeoMesh;
   firstmesh->NodeVec().Resize(4);
   TPZVec<REAL> coord(2);
   coord[1] = 0.;
   coord[0] = 0.;
   //nos geometricos
   firstmesh->NodeVec()[0].Initialize(coord,*firstmesh);
   coord[0] = 0.5;
   firstmesh->NodeVec()[1].Initialize(coord,*firstmesh);
   coord[0] = 1.;
   firstmesh->NodeVec()[2].Initialize(coord,*firstmesh);
   coord[1] = 0.5;
   firstmesh->NodeVec()[3].Initialize(coord,*firstmesh);
   TPZVec<int> nodeindexes(3);
   nodeindexes[0] = 0;
   nodeindexes[1] = 1;
   nodeindexes[2] = 2;
   //elementos geometricos
   TPZGeoEl1d *el1 = new TPZGeoEl1d(nodeindexes,1,*firstmesh);
   new TPZGeoEl1d(nodeindexes,1,*firstmesh);
	new TPZGeoEl1d(nodeindexes,1,*firstmesh);
   nodeindexes.Resize(2);
   nodeindexes[0] = 2;
   nodeindexes[1] = 3;
	new TPZGeoEl1d(nodeindexes,1,*firstmesh);//4 elementos
   //Arquivos de saidas
	ofstream outgm1("gm1info.dat");
	ofstream outgm2("gm2info.dat");
	ofstream outgm3("gm3info.dat");
   ofstream outcm1("cm1info.dat");
	ofstream outcm2("cm2info.dat");
   //montagem de conectividades entre elementos
   firstmesh->BuildConnectivity();
	//teste de divisao geometrica
 	firstmesh->Print(outgm1);
   outgm1.flush();
   TPZVec<TPZGeoEl *> sub(2);
   TPZAdmChunkVector<TPZGeoEl *> &elements = firstmesh->ElementVec();
   //fim teste divisao
   //malha computacional
   TPZCompMesh *secondmesh = new TPZCompMesh(firstmesh);
   //material
   int matindex = secondmesh->MaterialVec().AllocateNewElement();
   TPZFMatrix b(1,1,1.),c(1,1,0.);
   TPZMat1dLin *mat = new TPZMat1dLin(1);
   mat->SetMaterial(c,c,b,b);
   secondmesh->MaterialVec()[matindex] = mat;
   //CC : condicao de contorno
 	TPZGeoElBC bcel1(el1,2,-1);
   int bcindex = firstmesh->BCElementVec().AllocateNewElement();
   firstmesh->BCElementVec()[bcindex] = bcel1;
	TPZFMatrix val1(1,1,0.), val2(1,1,1.);
   TPZBndCond *bc = mat->CreateBC(-1,0,val1,val2);
   bcindex = secondmesh->BndCondVec().AllocateNewElement();
   secondmesh->BndCondVec()[bcindex] = bc;
   //constroe a malha computacional
   secondmesh->AutoBuild();
   //test ExpandConnected
//fim test ExpandConnected
	//Divisao aleatoria
   //numeros minimos e maximos de elementos que a malha nao devera ultrapassar apos agrupamento e divisao, respectivamente
   CycleRefinements(*secondmesh,20,50);
	//Divisao manual
/*   TPZManVector<int> subindex;
   secondmesh->Divide(0,subindex);
   int ind = subindex[0];
	secondmesh->Divide(1,subindex);
   secondmesh->Divide(ind,subindex);
   secondmesh->Divide(subindex[0],subindex);
   firstmesh->Print(outgm1);
   outgm1.flush();
   secondmesh->Print(outcm1);
   //test HigherDimensionConnected : Positivo 1d
   //test SmallConnect
   secondmesh->Coarsen(subindex,subindex[0]);   */
   //fim divisao manual
   secondmesh->InitializeBlock();
   secondmesh->ComputeConnectSequence();
   secondmesh->Print(outcm1);
   outcm1.flush();
   firstmesh->Print(outgm1);
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
//   TPZElementMatrix ek,ef;
//   secondmesh->ElementVec()[0]->CalcStiff(ek,ef);
//	ek.fMat->Print();
//   ef.fMat->Print();
   delete secondmesh;
   delete firstmesh;
   return 0;
}
//FIM MAIN FIM MAIN FIM MAIN FIM MAIN FIM MAIN FIM MAIN FIM MAIN FIM MAIN
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

/////////////////////////////////////////////
//TESTE
/*   TPZManVector<TPZCompEl *> sub = 0;
   secondmesh->ElementVec()[0]->Divide(sub);
   delete secondmesh->ElementVec()[0];
   secondmesh->ElementVec()[0] = 0;
   secondmesh->ElementVec().SetFree(0);
   secondmesh->InitializeBlock();   */

/*   secondmesh->ElementVec()[1]->Refine(1,sub);

   secondmesh->InitializeBlock();
   firstmesh->Print(outgeo);
   secondmesh->Print(outcomp);     */
/////////////////////////////////////////////

//	virtual void Refine(int index,TPZManVector<TPZCompEl *> sub);

/*inline void TPZCompEl::Refine(int index,TPZManVector<TPZCompEl *> sub) {
	cout << "TPZCompEl::Refine is called." << endl;
} */

/*
   /////////////////////////TESTE : Cedric
   //TPZVec<int> permute(secondmesh->Block().NBlocks());
   //int i = 0;
   //while (i < permute.NElements()) {permute[i] = 8-i; i++;}
   //secondmesh->Permute(permute);
   /////////////////////////

*/
   //test HigherDimensionConnected
/*   int id = 0,side,targetdimension;
   TPZCompElSide store[8];
while (id != -1) {//teste
   cout << "\nEntre id do elemento geometrico e seu lado side \n";
   cin >> id >> side;
   cout << "Entre targetdimension\n";
   cin >> targetdimension;
   TPZGeoEl1d *gel = ( TPZGeoEl1d *) firstmesh->ElementVec()[id];
   TPZGeoElSide gelside(gel,side);
   TPZStack<TPZCompElSide> elsidevec(store,8);
   gelside.HigherDimensionConnected(targetdimension,elsidevec);
   if (elsidevec.NElements()) {
	   TPZCompElSide el = elsidevec[0];
	   TPZGeoEl *geoel = el.Element()->Reference();
	   cout << "id of geoel de side de dimensao maior = " << geoel->Id() << endl;
   } else { cout << "Elemento nao achado\n\n";}
}//teste
   //fim test HigherDimensionConnected  */

/*
   //test SmallConnect
   int id = 0,side,targetlevel;
   TPZCompElSide store[10];
while (id != -1) {//teste
   cout << "\nEntre id do elemento geometrico e seu lado side \n";
   cin >> id >> side;
   cout << "Entre targetlevel de menor level a ser barrida\n";
   cin >> targetlevel;
   TPZGeoEl1d *gel = ( TPZGeoEl1d *) firstmesh->ElementVec()[id];
   TPZGeoElSide gelside(gel,side);
   TPZStack<TPZCompElSide> elsidevec(store,10);
   gelside.SmallConnect(targetlevel,elsidevec,0);
   for (int i=0;i<elsidevec.NElements();i++) {
	   TPZCompElSide el = elsidevec[i];
	   TPZGeoEl *geoel = el.Element()->Reference();
	   cout << "id of geoel de side de dimensao maior = " << geoel->Id() << endl;
   }
   if(!elsidevec.NElements()) { cout << "Elemento nao achado\n\n";}
}//fim test SmallConnect

*/
/*
/*
   TPZGeoEl *geo;
   int index = 0;
   while (index >= 0) {
   	cout << "Entre index geometrico\n";
      cin >> index;
      if (index < 0) continue;
      elements = firstmesh->ElementVec();//atualiza mesh
      if (index > elements.NElements()) {
      	cout <<  "Indice " << index << " de elemento nao existente" << "\n";
			cout <<  "* * * No division * * *" << "\n";
         cout <<  "Indice maximo " << elements.NElements() << "\n\n";
         continue;
      }
		geo = elements[index];
		geo->Divide(sub);
      cout << "Elemento de id " << geo->Id() << " Dividido\n";
      cout << "Elementos obtidos de ids : " << sub[0]->Id() << " , " << sub[1]->Id() << "\n";
   	firstmesh->Print(outgm1);
	   outgm1.flush();
   }
*/

/*
   int id = 0,side;
   TPZCompElSide store[10],expand[10];
while (id != -1) {//teste
   cout << "\nEntre id do elemento geometrico e seu lado side \n";
   cin >> id >> side;
   TPZGeoEl1d *gel = ( TPZGeoEl1d *) firstmesh->ElementVec()[id];
   TPZCompElSide compside(gel->Reference(),side);
   TPZStack<TPZCompElSide> elsidevec(store,10);
   compside.HigherLevelElementList(elsidevec,0,1);
   TPZGeoElSide gelside = compside.Reference();
	TPZStack<TPZCompElSide> expandvec(expand,10);
   compside.ExpandConnected(elsidevec,1);
   for (int j=0;j<elsidevec.NElements();j++) {
	   TPZCompElSide el = elsidevec[j];
	   TPZGeoEl *geoel = el.Element()->Reference();
	   cout << "id/side de elsidevec = " << geoel->Id() << "/" << el.Side() << endl;
   }
   cout << "\n\n";
   for (int i=0;i<expandvec.NElements();i++) {
	   TPZCompElSide el = expandvec[i];
      TPZGeoEl *geoel;
	   if(el.Element()) {
           geoel = el.Element()->Reference();
           cout << "id/side de expandvec = " << geoel->Id() << "/" << el.Side() << endl;
           TPZTransform t(el.Reference().Dimension());
           el.Reference().SideTransform2(gelside,t);
           t.Mult().Print("multiplication matrix");
           t.Sum().Print("sum matrix");
           cout <<  "\nnumero ";
           cin >> side;//para quebrar a saida
      }
   }
   if(!expandvec.NElements()) { cout << "expandvec vacio\n\n";}
}
*/
/*   elements[0]->Divide(sub);
   elements[1]->Divide(sub);
   elements[4]->Divide(sub);
   elements[5]->Divide(sub);
   elements[10]->Divide(sub);
   elements[11]->Divide(sub);
   elements[9]->Divide(sub);  */

/*   elements[0]->Divide(sub);
   elements[1]->Divide(sub);
   elements[4]->Divide(sub);
   elements[5]->Divide(sub);
   elements[3]->Divide(sub);
   elements[11]->Divide(sub);
   elements[7]->Divide(sub);
   elements[16]->Divide(sub); */
 /*  elements[0]->Divide(sub);
   elements[1]->Divide(sub);
   elements[4]->Divide(sub);
   elements[5]->Divide(sub);
   elements[10]->Divide(sub);
   elements[11]->Divide(sub);
   elements[8]->Divide(sub);//aqui é 9
   elements[2]->Divide(sub);
   elements[19]->Divide(sub);
   elements[20]->Divide(sub);
   elements[21]->Divide(sub); */
