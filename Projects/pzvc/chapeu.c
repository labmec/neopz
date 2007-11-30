

#include "pzfmatrix.h"
#include "pzskylmat.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzgnode.h"
#include "pzsolve.h"
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

#include "pzmat2dlin.h"
#include "pzanalysis.h"
#include "pzmetis.h"
#include "pzplaca.h"

#include <stdio.h>
#include <time.h>
//template<class T>
//class TPZVec;
//#define NOTDEBUG

TPZMaterial *LerMaterial(char *filename);
void PressaoHid(TPZVec<REAL> &x,TPZVec<REAL> &force);
void ReadMesh(TPZGeoMesh &gmesh, istream &file);
void ReadMaterial(TPZCompMesh &cmesh, istream &file);



int main() {

   //malha geometrica
   TPZGeoMesh *firstmesh = new TPZGeoMesh;
   firstmesh->SetName("Malha Geometrica : Chapeu");
   ifstream file("chapeuani.dat");
   ReadMesh(*firstmesh,file);

   //Arquivos de saida
   ofstream outgm1("outgm1.dat");
   ofstream outcm1("outcm1.dat");
   ofstream outcm2("outcm2.dat");
   firstmesh->Print(outgm1);
   outgm1.flush();

   //montagem de conectividades entre elementos
   firstmesh->BuildConnectivity();
   //malha computacional
   TPZCompMesh *secondmesh = new TPZCompMesh(firstmesh);
   secondmesh->SetName("Malha Computacional : Chapeu");


   //material
   ReadMaterial(*secondmesh,file);
   //ordem de interpolacao
   int ord;
   cout << "Entre ordem 1,2,3,4,5 : ";
   cin >> ord;
//   TPZCompEl::gOrder = ord;
   firstmesh.SetDefaultOrder(order);
   //constru��o malha computacional
   TPZVec<int> csub(0);
   TPZManVector<TPZGeoEl *> pv(4);
   int n1=1,level=0;
   cout << "\nDividir ate nivel ? ";
   int resp;
   cin >> resp;
   int nelc = firstmesh->ElementVec().NElements();
   int el;
   TPZGeoEl *cpel;
   for(el=0;el<firstmesh->ElementVec().NElements();el++) {
     cpel = firstmesh->ElementVec()[el];
     if(cpel && cpel->Level() < resp)
		cpel->Divide(pv);

   }
   //analysis
   secondmesh->AutoBuild();
   firstmesh->Print(outgm1);
   outgm1.flush();
   secondmesh->AdjustBoundaryElements();
   secondmesh->InitializeBlock();
   secondmesh->Print(outcm1);
   TPZAnalysis an(secondmesh,outcm1);
   int numeq = secondmesh->NEquations();
   secondmesh->Print(outcm1);
   outcm1.flush();
   TPZVec<int> skyline;
   secondmesh->Skyline(skyline);
   TPZSkylMatrix *stiff = new TPZSkylMatrix(numeq,skyline);
   an.SetMatrix(stiff);
   an.Solver().SetDirect(ECholesky);
   secondmesh->SetName("Malha Computacional :  Connects e Elementos");
   // Posprocessamento
   an.Run(outcm2);
   TPZVec<char *> scalnames(5);
   scalnames[0] = "Mn1";
   scalnames[1] = "Mn2";
   scalnames[2] = "Vn1";
   scalnames[3] = "Vn2";
   scalnames[4] = "Deslocz";
   TPZVec<char *> vecnames(0);
   char plotfile[] =  "chapeu.pos";
   char pltfile[] =  "chapeu.plt";
   an.DefineGraphMesh(2, scalnames, vecnames, plotfile);
   an.Print("FEM SOLUTION ",outcm1);
   an.PostProcess(3);
   an.DefineGraphMesh(2, scalnames, vecnames, pltfile);
   an.PostProcess(2);
   firstmesh->Print(outgm1);
   outgm1.flush();
   delete secondmesh;
   delete firstmesh;
   return 0;
}

// ****************** fim do programa principal ****************

// rotina para calculo de pressao hidrostatica parede
// calcula pressao hidrostatica aplicada no plano vertical xz

void PressaoHid(TPZVec<REAL> &x,TPZVec<REAL> &force){
	force[2] = -1.+x[1];
}

// ****************** fim da rotina PressaoHid *****************


void ReadMesh(TPZGeoMesh &mesh,istream &file) {

	int nnod,nelem,nbound;
	file >> nnod >> nelem >> nbound;
	mesh.NodeVec().Resize(nnod);

	TPZVec<REAL> coord(3);
	int in,nodindex;
	for(in=0; in<nnod; in++) {
		file >> nodindex >> coord[0] >> coord[1] >> coord[2];
		mesh.NodeVec()[in].Initialize(coord,mesh);
	}
	TPZVec<int> nodeindexes(3);
	int elindex, matindex;
	for(in=0; in<nelem; in++) {
		file >> elindex >> matindex >> nodeindexes[0] >> nodeindexes[1] >> nodeindexes[2];
		TPZGeoElT2d *elg0 = new TPZGeoElT2d(nodeindexes,matindex,mesh);
	}
	int boundindex,side;
	for(in=0; in<nbound; in++) {
		file >> elindex >> side >> boundindex;
	TPZGeoElBC(mesh.ElementVec()[elindex],side,boundindex,mesh);
	}
}

void ReadMaterial(TPZCompMesh &mesh, istream &file) {
	int nmat,nbound;
	file >> nmat >> nbound;
	TPZFMatrix naxes(3,3);
	REAL ni1,ni2,h,E1,E2,G12,G13,G23,f;
	TPZVec<REAL> xf(6);
	int matindex;
	file >> f  >>  h  >>
   		   E1  >>  E2 >>
           G12 >> G13 >> G23 >>
           ni1 >> ni2;
	int ind;
	TPZPlaca *mat = 0;
	for(ind = 0; ind<nmat; ind++) {
		file >> matindex;
		file >> naxes(0,0) >> naxes(0,1) >> naxes(0,2) >> naxes(1,0)
			>> naxes(1,1) >> naxes(1,2) >> naxes(2,0) >> naxes(2,1)
			>> naxes(2,2);
		file >> xf[0] >> xf[1] >> xf[2] >> xf[3] >> xf[4] >> xf[5];

		mat = new TPZPlaca(matindex,h,f,E1,E2,ni1,ni2,G12,G13,G23,naxes,xf);
		mesh.InsertMaterialObject(mat);
	}
	int i,j, boundtype,boundindex;
	TPZFMatrix val1(6,6,0.),val2(6,1,0.);
	for(ind=0; ind<nbound; ind++){
		file >> boundindex >> boundtype;
		for(i=0; i<6; i++) for(j=0; j<6; j++) file >> val1(i,j);
		for(i=0; i<6; i++) file >> val2(i,0);
		TPZBndCond *bc;
		bc = mat->CreateBC(boundindex,boundtype,val1,val2);
		mesh.InsertMaterialObject(bc);
	}

}
