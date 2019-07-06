/***************************************************************************
                          main.c  -  description
                             -------------------
    begin                : Thu Jun 15 2000
    copyright            : (C) 2000 by Edimar Cesar Rylo
    email                : ecrylo@fec.unicamp.br
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include <iostream>
#include <cstdlib>
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzgeoel.h"
#include "pzelgq2d.h"
#include "pzmat2dlin.h"
#include "pzelgq2dcyl.h"
#include "pzdxmesh.h"
#include "tpzmultcamada.h"
#include "tpzmatplacaiso.h"
#include "pzskylmat.h"
#include "pzanalysis.h"
#include "pzsolve.h"
#include "pzcylinsys.h"

void CycleRefinements(TPZCompMesh&cm, int numcycles, int minel, int maxel, ofstream &out);
int IsGroup(int ind, TPZCompEl *cel, TPZVec<int> &indgroup);
void LerMaterial(TPZCompMesh &cmesh,char *filename,TPZMultCamada * &mat);


int main() {

  const REAL pi = 3.14159265359;
	double coordstore[12][3] = {{1.,0.,0.},{1.,pi/2.,0.},{1.,pi,0.},{1.,3.*pi/2.,0.},{1.,0.,1.},{1.,pi/2.,1.},{1.,pi,1.},{1.,3.*pi/2.,1.},
	  {1.,0.,2.},{1.,pi/2.,2.},{1.,pi,2.},{1.,3.*pi/2.,2.}};
	// criar um objeto tipo malha geometrica
	TPZGeoMesh malha;
	TPZCylinsys *cyl = new TPZCylinsys;
	int sysindex = malha.CosysVec().AllocateNewElement();
	malha.CosysVec()[sysindex]=cyl;

	// criar quatro nos
	int i,j;
	TPZVec<REAL> coord(3,0.);
	  for(i=0; i<12; i++) {
	    // initializar as coordenadas do no em um vetor
	    for (j=0; j<3; j++) coord[j] = coordstore[i][j];
	    cyl->ToCart(&coord[0]);

	    // identificar um espaço no vetor onde podemos armazenar
	    // este vetor
	    int nodeindex = malha.NodeVec ().AllocateNewElement ();

	    // initializar os dados do nó
	    malha.NodeVec ()[i].Initialize (i,coord,malha);
	  }

	  // criar um elemento
	int el;
	TPZGeoEl *gel[8];
	for(el=0; el<4; el++) {
	  // initializar os indices dos nós
	  TPZVec<int> indices(4);
//	  for(i=0; i<4; i++){
	  if (el==3 || el==7){
	    indices[0] = el;
	    indices[1] = el-3;
	    indices[2] = el+1;
	    indices[3] = el+4;
	  }
	  else {
	  indices[0] = el;
	  indices[1] = el+1;
	  indices[2] = el+5;
	  indices[3] = el+4;
	  }

	  // O proprio construtor vai inserir o elemento na malha
	  gel[el] = new TPZGeoElQ2dCyl(el,indices,1,malha,sysindex);
	  //gel[el] = new TPZGeoElQ2d(el,indices,1,malha);
	}
 	ofstream output("output.dat");

	//	malha.Print(output);
	//output.flush();
	malha.BuildConnectivity ();
	//malha.Print(output);

	TPZCompMesh comp(&malha);

	// inserir os materiais
//	TPZMatPlacaIso *meumat = new TPZMatPlacaIso();
//	ifstream matarq("/home/pos/cesar/pzrepository/project/multiplaca/multiplaca3.dat");
	TPZMultCamada *meumat;
	LerMaterial(comp,"multplaca2.dat",meumat);
//	TPZCompEl::gOrder = 4;
  comp.SetDefaultOrder(order);

	comp.AutoBuild();
	comp.AdjustBoundaryElements();
	comp.InitializeBlock();

	comp.Print(output);
	TPZAnalysis an(&comp,output);
	int numeq = comp.NEquations();
	TPZVec<int> skyline;
	comp.Skyline(skyline);
	// TPZFMatrix *stiff =  new TPZFMatrix(numeq,numeq);
	TPZSkylMatrix *stiff = new TPZSkylMatrix(numeq,skyline);
	an.SetMatrix(stiff);

	an.Solver().SetDirect(ECholesky);
	comp.SetName("Malha Computacional :  Connects e Elementos");
   // Posprocessamento
	an.Run(output);
	//	an.Solution()->Zero();
	//	an.LoadSolution();
	TPZVec<char *> scalnames(4);
	scalnames[0] = "N1Scal";
	scalnames[1] = "N2Scal";
	scalnames[2] = "M1Scal";
	scalnames[3] = "M2Scal";

	TPZVec<char *> vecnames(5);
	vecnames[0] = "Displacement";
	vecnames[1] = "N1Vec";
	vecnames[2] = "N2Vec";
	vecnames[3] = "M1Vec";
	vecnames[4] = "M2Vec";
	char plotfile[] =  "placaPos.dx";
	an.DefineGraphMesh(2, scalnames, vecnames, plotfile);
	an.Print("FEM SOLUTION ",output);
	an.PostProcess(3);

	return 0;
}

void CycleRefinements(TPZCompMesh& cm, int numcycles, int minel, int maxel, ofstream &out){

	int numel = cm.NElements();
	if(numel > minel) minel = numel;

	int ic;
	for(ic=0; ic<numcycles; ic++) {
		numel = cm.NElements()-cm.ElementVec().NFreeElements();
		while(numel < maxel) {
			int elindex = rand()%cm.NElements();
			TPZCompEl *cel = cm.ElementVec()[elindex];
			if(cel) {
				TPZVec<int> subindex;
				cel->Divide(elindex,subindex);
				numel = cm.NElements()-cm.ElementVec().NFreeElements();
			}
		}
		//cm.Print(out);
		int i;
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
		//out.flush();
		//cm.Print(out);
		cm.CleanUpUnconnectedNodes();
	}
}

int IsGroup(int ind, TPZCompEl *cel, TPZVec<int> &indgroup) {
	TPZCompMesh *msh = cel->Mesh();
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

void LerMaterial(TPZCompMesh &cmesh,char *filename, TPZMultCamada * &matcamadas) {
   ifstream input(filename);
   TPZFMatrix naxes(3,3);
   REAL ni1,ni2,h,E1,E2,G12,G13,G23,f;
   REAL n00,n01,n02,n10,n11,n12,n20,n21,n22;
   TPZVec<REAL> xf(6), esp(1);
   int matindex, numcamadas,i,numbc,camadaref, materialtype;
   input >> matindex >> numcamadas >> camadaref >> numbc >> materialtype;
   matcamadas = new TPZMultCamada(matindex);
   cmesh.InsertMaterialObject(matcamadas);
   esp.Resize(numcamadas);
   for(i=0; i<numcamadas; i++) input >> esp[i];
   for(i=0; i< numcamadas; i++) {
     input >> f  >>  h  >>
       E1  >>  E2 >>
       G12 >> G13 >> G23 >>
       ni1 >> ni2;
     input >> n00 >> n01 >> n02 >> n10 >> n11 >> n12 >> n20 >> n21 >> n22;
     if(materialtype == 1) {
       xf.Resize(3*(numcamadas+1));
       int i;
       for(i=0; i< 3*(numcamadas+1); i++) input >> xf[i];
     } else {
       input >> xf[0] >> xf[1] >> xf[2] >> xf[3] >> xf[4] >> xf[5];
     }
     naxes(0,0) =  n00;    naxes(0,1) =  n01;    naxes(0,2) =  n02;
     naxes(1,0) =  n10;    naxes(1,1) =  n11;    naxes(1,2) =  n12;
     naxes(2,0) =  n20;    naxes(2,1) =  n21;    naxes(2,2) =  n22;
//     if(materialtype == 1) {
//       matcamadas->AddLayer(new TPZMatPlacaIso(-1,h,esp,f,E1,E2,ni1,ni2,G12,G13,G23,naxes,xf,camadaref,i));
//     } else {
       matcamadas->AddLayer(new TPZMatPlacaIso(-1,h,f,E1,E2,ni1,ni2,G12,G13,G23,naxes,xf));
//     }
   }
   int belindex, elside, bnum, btype, idf, jdf;
   TPZFMatrix val1(matcamadas->NStateVariables(),matcamadas->NStateVariables());
   TPZFMatrix val2(matcamadas->NStateVariables(),1);
   for(i=0; i<numbc; i++) {
     input >> belindex >> elside;
     input >> bnum >> btype;
     for(idf=0; idf< matcamadas->NStateVariables(); idf++) {
       for(jdf=0; jdf< matcamadas->NStateVariables(); jdf++) {
 	input >> val1(idf,jdf);
       }
     }
     for(idf=0; idf< matcamadas->NStateVariables(); idf++) {
       input >> val2(idf,0);
     }
     TPZGeoMesh *gmesh = cmesh.Reference();
     TPZGeoEl *gel = gmesh->ElementVec()[belindex];
     TPZGeoElBC(gel,elside,bnum,*gmesh);
     cmesh.InsertMaterialObject(matcamadas->CreateBC(bnum,btype,val1,val2));
     //   mat = matcamadas;
   }
}
