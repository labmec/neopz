// -*- c++ -*-
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzmat2dlin.h"
#include "pzbndcond.h"
#include "pzdxmesh.h"

#include <fstream>

using std::ifstream;
using std::ofstream;

void LerMalha(char *arquivo,TPZGeoMesh &mesh);

int main() {

	double coordstore[4][3] = {{0.,0.,0.},{1.,0.,0.},{1.,1.,0.},{0.,1.,0.}};
	// criar um objeto tipo malha geometrica
	TPZGeoMesh malha;

	// criar quatro nos
	int i,j;
	TPZVec<REAL> coord(3,0.);
	malha.NodeVec().Resize(4);
	for(i=0; i<4; i++) {
		// initializar as coordenadas do no em um vetor
		for (j=0; j<3; j++) coord[j] = coordstore[i][j];

		// identificar um espa�o no vetor onde podemos armazenar
		// este vetor
		//		int nodeindex = malha.NodeVec ().AllocateNewElement ();

		// initializar os dados do n�
		malha.NodeVec ()[i].Initialize (i,coord,malha);
	}

	// criar um elemento

	// initializar os indices dos n�s
	TPZVec<int> indices(4);
	for(i=0; i<4; i++) indices[i] = i;

	// O proprio construtor vai inserir o elemento na malha
        int index;
	TPZGeoEl *gel = malha.CreateGeoElement(EQuadrilateral,indices,1,index,1);


	malha.BuildConnectivity ();

	// Associar o lado de um elemento com uma condicao de contorno
	// Este objeto ira inserir-se automaticamente na malha
	TPZGeoElBC(gel,4,-1,malha);

	malha.Print();

	TPZGeoMesh malha2;
	LerMalha("quad_st_800_1R_16X.gri",malha2);
	
	ofstream output("output.dat");
	malha2.Print(output);

	// criar malha computacional
	TPZCompMesh comp(&malha2);

	// inserir os materiais
        TPZMat2dLin *meumat  = new TPZMat2dLin(1);
	TPZAutoPointer<TPZMaterial> mat(meumat);
	TPZFMatrix xk(1,1,1.),xc(1,2,0.),xf(1,1,1.);
	meumat->SetMaterial (xk,xc,xf);
	comp.InsertMaterialObject(mat);

	// inserir a condicao de contorno
	TPZFMatrix val1(1,1,0.),val2(1,1,0.);
	TPZAutoPointer<TPZMaterial> bnd = mat->CreateBC (mat,-4,0,val1,val2);
	comp.InsertMaterialObject(bnd);

	comp.AutoBuild();
	comp.InitializeBlock();

	comp.Print(output);

	TPZVec<char *> scalarnames(1),vecnames(0);
	scalarnames[0] = "state";
	TPZDXGraphMesh graph(&comp,2,mat,scalarnames,vecnames);
	ofstream *dxout = new ofstream("output.dx");
	graph.SetOutFile(*dxout);
	graph.SetResolution(0);
 
	graph.DrawMesh(1);
	graph.DrawSolution(0,0.);
	return 0;

}

void LerMalha(char *nome, TPZGeoMesh &grid) {
	ifstream infile(nome);

	int linestoskip;
	char buf[256];
	infile >> linestoskip;
	int i,j;
	for(i=0; i<linestoskip;i++) infile.getline(buf,255);
	infile.getline (buf,255);
	infile.getline (buf,255);
	int ntri,npoin,nbouf,nquad,nsidif;
	infile >> ntri >> npoin >> nbouf >> nquad >> nsidif;
	infile.getline (buf,255);
	infile.getline(buf,255);

	grid.NodeVec ().Resize(npoin+1);
	TPZVec<int> nodeindices(4);
	int mat, elid;
	for(i=0;i<nquad;i++) {
		infile >> elid;
		for(j=0; j<4;j++) infile >> nodeindices[j];
		infile >> mat;
                int index;
		grid.CreateGeoElement(EQuadrilateral,nodeindices,mat,index,1);
	}
	infile.getline(buf,255);
	infile.getline(buf,255);

	int nodeid,dum;
	char c;
	TPZVec<REAL> coord(3,0.);
	for(i=0; i<npoin; i++) {
		infile >> nodeid >> coord[0] >> coord[1] >> c >> dum;
		grid.NodeVec ()[nodeid].Initialize (nodeid,coord,grid);
	}
	infile.getline (buf,255);
	infile.getline (buf,255);

	TPZVec<int> sideid(2,0);
	for(i=0; i<nbouf; i++) {
		infile >> sideid[0] >> sideid[1] >> elid >> dum >> mat;
		TPZGeoEl *el = grid.ElementVec ()[elid-1];
		int side = el->WhichSide (sideid);
		TPZGeoElBC(el,side,-mat,grid);
	}
	grid.BuildConnectivity();

	return;
}
