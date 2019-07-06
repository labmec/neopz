// viga3d.cpp : Defines the entry point for the console application.
//

#include "pzreal.h"
#include "pzvec.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzgeoel.h"
#include "pzelgc3d.h"
#include "pzmathyperelastic.h"
#include "pznonlinanalysis.h"
#include "pzdxmesh.h"
#include "pzskylmat.h"
#include "pzsolve.h"

#include <iostream.h>
#include <fstream.h>

int IsGroup(int ind, TPZCompEl *cel, TPZVec<int> &indgroup);

static REAL angle = 0.;
static REAL pi = 3.141592654;
void Forcing(TPZVec<REAL> &x, TPZVec<REAL> &disp);

int main(int argc, char* argv[]){

	cout << "Exemplo de aplicação de modelagem Tri-dimensional através do DX!\n";

	double coordstore[4][2] = {{0.,0.},{1.,0.},{1.,1.},{0.,1.}};
	int numel = 10;
	TPZVec<REAL> coord(3,0.);

	// criar um objeto tipo malha geometrica
	TPZGeoMesh malha;

	// criar nos
	int i,j;
	for(i=0; i<4*(numel+1); i++) {
		// initializar as coordenadas do no em um vetor
		for (j=0; j<2; j++) {
			coord[j] = coordstore[i%4][j];
		}
		coord[2] = i/4;
		// identificar um espaço no vetor onde podemos armazenar
		// este vetor
		int nodeindex = malha.NodeVec().AllocateNewElement();
		// initializar os dados do nó
		malha.NodeVec()[i].Initialize(i,coord,malha);
	}

	// criar um elemento
	int el;
	TPZGeoEl *gel[10];
	TPZVec<int> indices(8);

	for(el=0; el<numel; el++) {
		// inicializar os indices dos nós
		for(i=0; i<8; i++) indices[i] =4*el+i;
		// O proprio construtor vai inserir o elemento na malha
		gel[el] = new TPZGeoElC3d(el,indices,1,malha);
	}
/*	TPZVec<TPZGeoEl *> sub;
	gel[1]->Divide(sub);*/

	ofstream output("output3d.dat");

	malha.Print(output);
	output.flush();
	malha.BuildConnectivity();
		TPZVec<TPZGeoEl *> sub;
	gel[0]->Divide(sub);
	gel[9]->Divide(sub);
	// TPZMaterial::gBigNumber = 1.e6;

	TPZGeoElBC t3(gel[0],20,-1,malha);
/*	TPZGeoElBC(gel[0],0,-3,malha);
	TPZGeoElBC(gel[0],1,-4,malha);
	TPZGeoElBC(gel[0],2,-5,malha);
*/
	TPZGeoElBC t4(gel[numel-1],25,-2,malha);
	malha.Print(output);

	TPZCompMesh comp(&malha);

	// inserir os materiais
	TPZMaterial *meumat = new TPZMatHyperElastic(1,1.e5,0.25);
	comp.InsertMaterialObject(meumat);

	// inserir a condicao de contorno
	TPZFMatrix val1(3,3,0.),val2(3,1,0.);

	TPZMaterial *bnd = meumat->CreateBC (-1,0,val1,val2);
	comp.InsertMaterialObject(bnd);
	bnd = meumat->CreateBC (-2,0,val1,val2);
	bnd->SetForcingFunction(Forcing);
	comp.InsertMaterialObject(bnd);
/*	val2.Zero();
	bnd = meumat->CreateBC (-3,0,val1,val2);
	comp.InsertMaterialObject(bnd);
	val1(1,1) = 1.e10;
	val1(2,2) = 1.e10;
	bnd = meumat->CreateBC (-4,2,val1,val2);
	comp.InsertMaterialObject(bnd);
	val1.Zero();
	val1(2,2) = 1.e10;
	bnd = meumat->CreateBC (-5,2,val1,val2);
	comp.InsertMaterialObject(bnd);
*/
//	TPZCompEl::gOrder = 2;
  cmesh.SetDefaultOrder(2);
	comp.AutoBuild();

	angle = pi/4.;
	TPZNonLinearAnalysis an(&comp,output);
	int numeq = comp.NEquations();
	TPZVec<int> skyline;
	comp.Skyline(skyline);
	//TPZFMatrix *stiff = new TPZFMatrix(numeq,numeq);
	TPZSkylMatrix *stiff = new TPZSkylMatrix(numeq,skyline	);
	an.SetMatrix(stiff);
	an.Solution()->Zero();
	an.Solver().SetDirect(ELDLt);//ELDLt , ECholesky

	TPZDXGraphMesh graph(&comp,3,meumat);
	ofstream *dxout = new ofstream("output3d4.dx");
	graph.SetOutFile(*dxout);
	graph.SetResolution(2);

	graph.DrawMesh(2);
    TPZVec<char *> scalnames(1),vecnames(2);
    scalnames[0] = "VonMises";
    vecnames[0] = "state";
	vecnames[1] = "displacement";
	graph.DrawSolution(0,0.,scalnames,vecnames);

	int numiter =2;
	REAL tol = 1.e-5;
	an.IterativeProcess(output,tol,numiter);
	graph.DrawSolution(1,angle,scalnames,vecnames);
	//angle += pi/4.;
	//an.IterativeProcess(output,tol,numiter);
	//graph.DrawSolution(2,angle,scalnames,vecnames);
	/*
	angle += pi/4.;
	an.IterativeProcess(output,tol,numiter);
	graph.DrawSolution(3,angle,scalnames,vecnames);
	angle += pi/4.;
	an.IterativeProcess(output,tol,numiter);
	graph.DrawSolution(3,angle,scalnames,vecnames);
	angle += pi/4.;
	an.IterativeProcess(output,tol,numiter);
	graph.DrawSolution(3,angle,scalnames,vecnames);
	*/

	//char plotfile[] =  "output3d2.dx";
	//an.DefineGraphMesh(3, scalnames, vecnames, plotfile);

	//comp.Print(output);

	//graph.DrawSolution(0,0.,scalnames,vecnames);


	return 0;
}


void Forcing(TPZVec<REAL> &x, TPZVec<REAL> &disp){
	disp[0] = -(x[1]-0.5)*sin(angle)+(x[0]-0.5)*cos(angle)-(x[0]-0.5);
	disp[1] = (x[1]-0.5)*cos(angle)+(x[0]-0.5)*sin(angle)-(x[1]-0.5);
	disp[2] = 0.;
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



