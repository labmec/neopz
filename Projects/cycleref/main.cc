//#include "pzmetis.h"
//#include "pztrnsform.h"
#include "pzgeoel.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzfmatrix.h"
#include "pzskylmat.h"
#include "pzgnode.h"
#include "pzanalysis.h"
#include "pzcompel.h"
#include "pzbndcond.h"
#include "pzstack.h"
#include "pzvec.h"
#include "pzsolve.h"
#include "pzelg1d.h"
#include "pzelgq2d.h"
#include "pzelgt2d.h"
#include "pzelct2d.h"
#include "pzelgc3d.h"
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
#include "pzmattest3d.h"

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <iostream.h>

//#define NOTDEBUG

void Divide(TPZCompMesh *compmesh);
void LerMalha(char *malha,TPZGeoMesh *geomesh,TPZCompMesh *& compmesh);
void CycleRefinements(TPZCompMesh&cm, int numcycles, int minel, int maxel, ofstream &out);
int IsGroup(int ind, TPZCompEl *cel, TPZVec<int> &indgroup);

int main() {

	// criar um objeto tipo malha geometrica
	TPZGeoMesh malha;
   //Arquivos de saida
	ofstream outgm("outgm.txt");
   ofstream outcm("outcm.txt");
   TPZCompMesh *particao;
	// criar elementos
   LerMalha("malha.txt",&malha,particao);
   //montagem de conectividades entre elementos
	malha.BuildConnectivity();
   //ordem de interpolacao
   int ord;
   cout << "Enter order 1,2,3,4,5,... : \n";
   cin >> ord;
   TPZCompEl::gOrder = ord;
   //construção malha computacional
	particao->AutoBuild();
   //analysis
	particao->InitializeBlock();
   //malha antes do refinamento
   malha.Print(outgm);
   outgm.flush();
   particao->Print(outcm);
   outcm.flush();
   //divisão manual
   int nao;
   cout << "\nRemove manual (0/1)?";
   cin >> nao;
   int maxel,sim;
   ofstream output("output.dat");
   if(nao) Divide(particao);
   //teste de restrições
   else {
      cout << "\nEntre maxel : ";
      cin >> maxel;
      CycleRefinements(*particao,50,2,maxel,output);//maxel=150
   }
   //malha apos refinamento
	particao->Print(output);
   malha.Print(outgm);
   outgm.flush();
   particao->Print(outcm);
   outcm.flush();
   cout << "\nEnd & \n";
   cin >> sim;
   delete particao;
	return 0;
}

void CycleRefinements(TPZCompMesh& cm, int numcycles, int minel, int maxel, ofstream &/*out*/){

	int numel = cm.NElements();
	if(numel > minel) minel = numel;
   TPZAdmChunkVector<TPZCompEl *> &elemvec = cm.ElementVec();
	int ic;
   randomize();
   clock_t start,end,begin;
   begin = clock();
	for(ic=0; ic<numcycles; ic++) {
      start = clock();
      cout << "\nProcessed cycle number : " << (ic-1) << endl;
		numel = cm.NElements() - cm.ElementVec().NFreeElements();
		while(numel < maxel) {
			int elindex = rand()%cm.NElements();
			TPZCompEl *cel = elemvec[elindex];//cm.ElementVec()[elindex];
			if(cel) {
				TPZVec<int> subindex;
				cel->Divide(elindex,subindex);
				numel = cm.NElements() - cm.ElementVec().NFreeElements();
			}
		}
      end = clock();
      cout << "\nEnd divide cycle\n";
      cout << ((end - start)/CLK_TCK) << " segundos" << endl;
      start = end;
		//cm.Print(out);
		//int i;
		while(numel > minel) {
			int elindex = rand()%cm.NElements();
			TPZCompEl *cel = cm.ElementVec()[elindex];
			if(cel) {
				TPZVec<int> subindex;
				if(IsGroup(elindex,cel,subindex)){
				  //out << "Coarsening "; for(i=0;i<4;i++) out << subindex[i] << "/" << cm.ElementVec()[subindex[i]]->Reference()->Id() << " ";
				  //out << endl;
					cm.Coarsen(subindex,elindex);
					//out << "Created " << elindex << "/" << cm.ElementVec()[elindex]->Reference()->Id() << endl;
					//out.flush();
					cm.CleanUpUnconnectedNodes();
					numel = cm.NElements()-cm.ElementVec().NFreeElements();
				}
			}

		}
      end = clock();
      cout << "\nEnd coarsen cycle\n";
      cout << ((end - start)/CLK_TCK) << " segundos" << endl;
		//out.flush();
		//cm.Print(out);
		cm.CleanUpUnconnectedNodes();
      end = clock();
      cout << "\nEnd cycle : total time\n";
      cout << ((end - begin)/CLK_TCK) << " segundos" << endl;
	}
}

int IsGroup(int /*ind*/, TPZCompEl *cel, TPZVec<int> &indgroup) {
	//TPZCompMesh *msh = cel->Mesh();
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

void LerMalha(char *malha,TPZGeoMesh *geomesh,TPZCompMesh *&compmesh) {

	ifstream grid(malha);
   TPZFMatrix xk(1,1,0.),xb(1,1,0.),xc(1,1,0.),xf(1,1,0.);
   TPZMaterialTest3D *mat1d,*mat2d,*mat3d;//TPZMat1dLin *mat1d;//TPZMat2dLin *mat2d;
   int nnode,nel,ncorners;
   TPZVec<REAL> coord(3);
   grid >> nnode >> nel;
   geomesh->NodeVec().Resize(nnode);
   for(int inode=0;inode<nnode;inode++) {
      grid >> coord[0];
      grid >> coord[1];
      grid >> coord[2];
		// identificar um espaço no vetor onde podemos armazenar este ponto
/*int nodeindex = */geomesh->NodeVec().AllocateNewElement();
      geomesh->NodeVec()[inode].Initialize(inode,coord,*geomesh);
   }
   compmesh = new TPZCompMesh(geomesh);
   for(int el=0;el<nel;el++)  {
      int mat,ntype;
      grid >> ntype >> mat;
      ncorners = ntype;
      if(ntype == 7) ncorners = 4;
      TPZVec<int> nodes(ncorners);
      for(int incid=0;incid<ncorners;incid++) grid >> nodes[incid];
      switch(ntype) {//tipo de elemento
         case 2://unidimensional ; elg1d =
            new TPZGeoEl1d(nodes,mat,*geomesh);
            mat1d = new TPZMaterialTest3D(mat);//mat1d = new TPZMat1dLin(mat);
            compmesh->InsertMaterialObject(mat1d);
            mat1d->SetMaterial(xk);//mat1d->SetMaterial(xk,xb,xc,xf);
            break;
         case 3://triângulo ; elgt2d =
         	new TPZGeoElT2d(nodes,mat,*geomesh);
            mat2d = new TPZMaterialTest3D(mat);//mat2d = new TPZMat2dLin(mat);
            compmesh->InsertMaterialObject(mat2d);//mat2d->SetMaterial(xk,xc,xf);
            mat2d->SetMaterial(xk);
            break;
         case 4://quadrilátero ; elgq2d =
            new TPZGeoElQ2d(nodes,mat,*geomesh);
            mat2d = new TPZMaterialTest3D(mat);//mat2d = new TPZMat2dLin(mat);
            compmesh->InsertMaterialObject(mat2d);//mat2d->SetMaterial(xk,xc,xf);
            mat2d->SetMaterial(xk);
            break;
         case 7://tetraedro ; elgt3d =
            new TPZGeoElT3d(nodes,mat,*geomesh);
            mat3d = new TPZMaterialTest3D(mat);
            compmesh->InsertMaterialObject(mat3d);
            mat3d->SetMaterial(xk);
            break;
         case 5://pirâmide ; elgpi3d =
            new TPZGeoElPi3d(nodes,mat,*geomesh);
            mat3d = new TPZMaterialTest3D(mat);
				compmesh->InsertMaterialObject(mat3d);
            mat3d->SetMaterial(xk);
            break;
         case 6://pirâmide ; elgpi3d =
            new TPZGeoElPr3d(nodes,mat,*geomesh);
            mat3d = new TPZMaterialTest3D(mat);
				compmesh->InsertMaterialObject(mat3d);
            mat3d->SetMaterial(xk);
            break;
         case 8://cubo ; elgc3d =
            new TPZGeoElC3d(nodes,mat,*geomesh);
            mat3d = new TPZMaterialTest3D(mat);
				compmesh->InsertMaterialObject(mat3d);
            mat3d->SetMaterial(xk);
	         break;
         default:
         	for(int i=0;i<300;i++)
            	cout << "\nmain::LerMalha -> Elemento nao conhecido\n";
            cout << "\nChao\n";
            exit(1);
            delete mat1d;
            delete mat2d;
            delete mat3d;
      }
   }
}
void Divide(TPZCompMesh *compmesh) {

   TPZVec<int> csub;
   int n1=1,remove;
   while(n1) {
	   cout << "Id do elemento geometrico remover ? : ";
      cin >> n1;
      if(n1 < 0) break;
      cout << "\n<1:Divide> : <2:Agrupa> ";
      cin >> remove;
      int nelc = compmesh->ElementVec().NElements();
      int el;
      TPZCompEl *cpel;
      for(el=0;el<nelc;el++) {
         cpel = compmesh->ElementVec()[el];
         if(cpel && cpel->Reference()->Id() == n1) break;
         if(remove==1) compmesh->Divide(el,csub,0);
         if(remove==2) {
            if(IsGroup(el,cpel,csub))
               compmesh->Coarsen(csub,el);
         }
         compmesh->CleanUpUnconnectedNodes();
      }
      n1 = 1;
   }
}
