/**
 * @file
 * @brief Implements the build of a computational mesh as tutorial example of the interpolation NeoPZ module
 */
#include <pzvec.h>
#include <pzgmesh.h>
#include <pzcompel.h>
#include <pzgeoel.h>
#include <pzquad.h>
#include <pzmat2dlin.h>
#include <TPZGeoElement.h>
#include <pzskylstrmatrix.h>
#include <pzcmesh.h>
#include "pzfilebuffer.h"
#include "pzmaterialid.h"
#include "pzmeshid.h"
#include "pzbfilestream.h"
#include "pzcmesh.h"
#include <TPZVTKGeoMesh.h>
#include "pzelast3d.h"
#include "pzstepsolver.h"
#include "pzanalysis.h"
#include "tpzautopointer.h"
#include "pzmaterial.h"
#include "pzbndcond.h"

#include <iostream>
#include <fstream>
using namespace std;

// nx = number of nodes in x direction
// ny = number of nodes in y direction
TPZGeoMesh * GetMesh(int nx,int ny);
void SimpleMesh(TPZGeoMesh *malha);
void InsertElasticity(TPZCompMesh *cmesh);


int main()
{
	TPZGeoMesh *gMesh = new TPZGeoMesh;
	SimpleMesh(gMesh);
	TPZCompMesh *cmesh = new TPZCompMesh(gMesh);
	cmesh->SetDefaultOrder(3);
	cmesh->SetDimModel(3);
	InsertElasticity(cmesh);
	cmesh->AutoBuild();
	TPZSkylineStructMatrix skylstruct(cmesh);
	TPZStepSolver<REAL> step;
	step.SetDirect(ECholesky);
	TPZAnalysis an(cmesh);
	an.SetStructuralMatrix(skylstruct);
	an.SetSolver(step);
	an.Run();	
	
	return 0;
}


void SimpleMesh(TPZGeoMesh *malha)
{
    
	double coordstore[][3] = {{0.,0.,0.},{1.,0.,0.},{1.,1.,0.},{0.,1.,0.}
        , {0.,0.,1.},{1.,0.,1.},{1.,1.,1.},{0.,1.,1.}};
    
	// criar quatro nos
	int i,j;
	TPZVec<REAL> coord(3,0.);
	malha->NodeVec().Resize(8);
	for(i=0; i<8; i++) {
		// initializar as coordenadas do no em um vetor
		for (j=0; j<3; j++) coord[j] = coordstore[i][j];
        
		// identificar um espaco no vetor onde podemos armazenar
		// este vetor
		//		int nodeindex = malha.NodeVec ().AllocateNewElement ();
        
		// initializar os dados do no
		malha->NodeVec ()[i].Initialize (i,coord,*malha);
	}
    
	// criar um elemento
    
	// initializar os indices dos n�s
	TPZVec<int> indices(8);
	for(i=0; i<8; i++) indices[i] = i;
    
	// O proprio construtor vai inserir o elemento na malha
    int index;
	malha->CreateGeoElement(ECube,indices,1,index,0);

	malha->BuildConnectivity ();
	
	/*
	int iref;
	TPZManVector<TPZGeoEl *> subs;
	for (iref=0; iref<nrefloop; iref++) 
	{
		int nel = malha.NElements();
		int iel;
		for (iel=0; iel<nel; iel++) {
			malha.ElementVec()[iel]->Divide(subs);
		}
		std::cout << "Refinement step " << iref << " Number of elements " << malha.NElements() << std::endl;
		std::cout.flush();
	}
	 */
	int el, dirID = -1, numelements = malha->NElements();
	TPZManVector <int> TopolCubo(8);
	
	for (el=0; el<numelements; el++)
	{
		TPZGeoEl *cubo = malha->ElementVec()[el];
		for (int i=0; i<8; i++)
		{
			TopolCubo[i] = cubo->NodeIndex(i); 
		}
		
		// Colocando as condicoes de contorno
		TPZManVector <TPZGeoNode,4> Nodefinder(8);
		TPZManVector <REAL,3> nodecoord(3);
		// na face x = 1
		TPZVec<int> ncoordzVec(0); int sizeOfVec = 0;
		for (int i = 0; i < 8; i++) 
		{
			Nodefinder[i] = malha->NodeVec()[TopolCubo[i]];
			Nodefinder[i].GetCoordinates(nodecoord);
			if (nodecoord[2] == 0.)
			{
				sizeOfVec++;
				ncoordzVec.Resize(sizeOfVec);
				ncoordzVec[sizeOfVec-1] = TopolCubo[i];
			}
		}
		if(ncoordzVec.NElements() == 4)
		{
			int lado = cubo->WhichSide(ncoordzVec);
			TPZGeoElSide cuboside(cubo, lado);
			TPZGeoElBC(cuboside,dirID);		
		}
		
		ncoordzVec.Resize(0);
		sizeOfVec = 0;
	}
	
	
	ofstream tst("cubo.vtk");
	TPZVTKGeoMesh::PrintGMeshVTK(malha, tst, true);
}

void InsertElasticity(TPZCompMesh *mesh)
{
	mesh->SetDimModel(3);
	int nummat = 1, dirichlet = 0, neumann = 1, mixed = 2;
	int dir = -1;
	TPZManVector<REAL> force(3,0.);
	force[2] = 5.;
	REAL Ela = 1000, poisson = 0.2; 
    
	TPZElasticity3D *elast = new TPZElasticity3D(nummat, Ela, poisson, force);

	TPZAutoPointer<TPZMaterial> elastauto(elast);
	mesh->InsertMaterialObject(elastauto);
	
	// Neumann em x = 1;
	TPZFMatrix<REAL> val1(3,3,0.),val2(3,1,0.);
	TPZBndCond *bc = elast->CreateBC(elastauto, dir, dirichlet, val1, val2);
	TPZAutoPointer<TPZMaterial> bcauto(bc);
	mesh->InsertMaterialObject(bcauto);
	
	
	/*
	// Dirichlet em 1 -1 -1 yz;
	val1(0,0) = 0.;
	val1(1,1) = 1.;
	val1(2,2) = 1.;
	TPZBndCond *bc2 = viscoelast->CreateBC(viscoelastauto, dir2, mixed, val1, val2);
	TPZAutoPointer<TPZMaterial> bcauto2(bc2);
	mesh->InsertMaterialObject(bcauto2);
	
	// Dirichlet em 1 1 -1 z;
	val1(0,0) = 0.;
	val1(1,1) = 0.;
	val1(2,2) = 1.;
	TPZBndCond *bc3 = viscoelast->CreateBC(viscoelastauto, dir3, mixed, val1, val2);
	TPZAutoPointer<TPZMaterial> bcauto3(bc3);
	mesh->InsertMaterialObject(bcauto3);
	 */
}


TPZGeoMesh *GetMesh (int nx,int ny){
	int i,j;
	int id, index;
	
	//Let's try with an unitary domain
	REAL lx = 1.;
	REAL ly = 1.;
	
	//Creates the geometric mesh... The nodes and elements
	//will be inserted into mesh object during initilize process
	TPZGeoMesh *gmesh = new TPZGeoMesh();
	
	//Auxiliar vector to store a coordinate
	TPZVec <REAL> coord (3,0.);
	
	//Nodes initialization
	for(i = 0; i < nx; i++){
		for(j = 0; j < ny; j++){
			id = i*ny + j;
			coord[0] = (i)*lx/(nx - 1);
			coord[1] = (j)*ly/(ny - 1);
			//using the same coordinate x for z
			coord[2] = 0.;
			//cout << coord << endl;
			//Get the index in the mesh nodes vector for the new node
			index = gmesh->NodeVec().AllocateNewElement();
			//Set the value of the node in the mesh nodes vector
			gmesh->NodeVec()[index] = TPZGeoNode(id,coord,*gmesh);
		}
	}

	
	//Creates a vector of pointers to geometric elements
	TPZGeoEl * elvec[(const int)((nx-1)*(ny-1))];
	
	//Auxiliar vector to store a element connectivities
	TPZVec <int> connect(4,0);
	
	//Element connectivities
	for(i = 0; i < (nx - 1); i++){
		for(j = 0; j < (ny - 1); j++){
			index = (i)*(ny - 1)+ (j);
			connect[0] = (i)*ny + (j);
			connect[1] = connect[0]+(ny);
			connect[2] = connect[1]+1;
			connect[3] = connect[0]+1;
			//cout << connect << endl;
			//Allocates and define the geometric element
			elvec[index] = gmesh->CreateGeoElement(EQuadrilateral,connect,1,id);
		}
	}
	//Generate neighborhod information
	gmesh->BuildConnectivity();
	ofstream tst("teste.vtk");
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, tst);
	return gmesh;
}