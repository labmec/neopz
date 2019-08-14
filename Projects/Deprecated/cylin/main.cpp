
#include <iostream.h>
#include <stdlib.h>
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzgeoel.h"
#include "pzelgq2d.h"
#include "pzmat2dlin.h"
#include "pzelgq2dcyl.h"
#include "pzdxmesh.h"

void CycleRefinements(TPZCompMesh&cm, int numcycles, int minel, int maxel, ofstream &out);
int IsGroup(int ind, TPZCompEl *cel, TPZVec<int> &indgroup);

int main() {

  const REAL pi = 3.14159265359;
	double coordstore[6][3] = {{1.,0.,0.},{1.,2.*pi/3.,0.},{1.,4.*pi/3.,0.},
	  {1.,0.,1.},{1.,2.*pi/3.,1.},{1.,4.*pi/3.,1.}};
	// criar um objeto tipo malha geometrica
	TPZGeoMesh malha;
	TPZCylinsys *cyl = new TPZCylinsys;
	int sysindex = malha.CosysVec().AllocateNewElement();
	malha.CosysVec()[sysindex]=cyl;

	// criar quatro nos
	int i,j;
	TPZVec<REAL> coord(3,0.); 
	  for(i=0; i<6; i++) {
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
	TPZGeoEl *gel[3];
	for(el=0; el<3; el++) {
	  
	  // initializar os indices dos nós
	  TPZVec<int> indices(4);
	  for(i=0; i<4; i++) indices[i] = (el+i%2)%3+(i/2)*3;
	  int tmp = indices[3];
	  indices[3]=indices[2];
	  indices[2]=tmp;

	  // O proprio construtor vai inserir o elemento na malha
	  gel[el] = new TPZGeoElQ2dCyl(el,indices,1,malha,sysindex);
	}
 	ofstream output("output.dat");

	malha.Print(output);
	output.flush();
	malha.BuildConnectivity ();
	malha.Print(output);

	TPZCompMesh comp(&malha);

	// inserir os materiais
	TPZMat2dLin *meumat = new TPZMat2dLin(1);
	TPZFMatrix xk(1,1,1.),xc(1,2,0.),xf(1,1,1.);
	meumat->SetMaterial (xk,xc,xf);
	comp.InsertMaterialObject(meumat);

	// inserir a condicao de contorno
	TPZFMatrix val1(1,1,0.),val2(1,1,0.);
	TPZMaterial *bnd = meumat->CreateBC (-4,0,val1,val2);
	comp.InsertMaterialObject(bnd);

	comp.AutoBuild();
	comp.InitializeBlock();

	comp.Print(output);

	TPZDXGraphMesh graph(&comp,2,meumat);
	ofstream *dxout = new ofstream("output.dx");
	graph.SetOutFile(*dxout);
	graph.SetResolution(5);
 
	graph.DrawMesh(1);
	TPZVec<char *> scalarnames(1),vecnames(0);
	scalarnames[0] = "state";
	graph.DrawSolution(0,0.,scalarnames,vecnames);
	/*
	TPZVec<REAL> coord1(3,0.);
	TPZVec<REAL> normal(3,0.);
	TPZFMatrix jacobian(2,2,0.);
	TPZFMatrix jacinv(2,2,0.);
	TPZFMatrix jac1d(1,1,0.);
	TPZFMatrix axes (3,3,0.);
	cyl->GetAxes(axes);
	REAL detjac;
	int k=0;
	for (k=0;k<3;k++){
	  gel[k]->Jacobian(coord1,jacobian,axes,detjac,jacinv);
	  for (i=0;i<2;i++)
	    for (j=0;j<2;j++){
	      output << "Jacobiano\t" << i << " " << j << "\t" << jacobian(i,j) << "\n";
	    }
	  TPZVec<REAL> loc(3,0.),result(3,0.);
	  gel[k]->X(loc,result);
	  output << "X\t" << result[0] <<"\t "<<result[1] << "\t" << result[2] <<"\n";
	}
	TPZVec<REAL> zero(3,0.);
	for (k=0;k<3;k++){
	  cyl->GetAxes(axes);
	  gel[k]->NormalVector(3,zero,normal,axes,jac1d);
	  for (i=0;i<3;i++){
	      output << "Normal\t" << i <<"\t" << normal[i] << "\n";
	  }
	}
	*/
	// criar malha computacional
// 	TPZCompMesh comp(&malha);

// 	// inserir os materiais
// 	TPZMat2dLin *meumat = new TPZMat2dLin(1);
// 	TPZFMatrix xk(1,1,1.),xc(1,2,0.),xf(1,1,1.);
// 	meumat->SetMaterial (xk,xc,xf);
// 	comp.InsertMaterialObject(meumat);

// 	// inserir a condicao de contorno
// 	TPZFMatrix val1(1,1,0.),val2(1,1,0.);
// 	TPZMaterial *bnd = meumat->CreateBC (-4,0,val1,val2);
// 	comp.InsertMaterialObject(bnd);

// 	comp.AutoBuild();
// 	comp.InitializeBlock();

// 	CycleRefinements(comp,50,2,150,output);
// 	comp.Print(output);

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
