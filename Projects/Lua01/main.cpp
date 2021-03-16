/**
 * @file
 * @brief Implements the build of a computational mesh as tutorial example of the interpolation NeoPZ module
 */

#include "pzcmesh.h"
#include <TPZVTKGeoMesh.h>
#include "pzelast3d.h"
#include "pzstepsolver.h"
#include "pzanalysis.h"
#include "tpzautopointer.h"
#include "TPZMaterial.h"
#include "pzbndcond.h"
#include <pzvec.h>
#include <pzgmesh.h>
#include <pzcompel.h>
#include <pzgeoel.h>
#include <pzquad.h>
#include <pzmat2dlin.h>
#include <TPZGeoElement.h>
#include <pzskylstrmatrix.h>
#include <pzcmesh.h>

#include "pzgeoelbc.h"

#include <pzelast3d.h>
#include <pzplaca.h>
#include <pzvtkmesh.h>
//#include <pzlog.h>


#include <iostream>
#include <fstream>

using namespace std;

// nx = number of nodes in x direction
// ny = number of nodes in y direction
TPZGeoMesh * GetMesh(int nx,int ny);
void InsertElasticity(TPZCompMesh *cmesh);

int main(){
	
	
	//  RegisterMeshClasses();
	//  RegisterMatrixClasses();
	//  RegisterMaterialClasses();
	int nx = 5,ny = 5;
	int order = 3;
	TPZCompEl::SetgOrder(order);
	
	//Creates the geometric mesh
	TPZGeoMesh *mesh = GetMesh(nx,ny);
	mesh->SetName("testing a space");
	ofstream out("all.dat");

	// Print geometrical mesh info
	ofstream geomeshout("geomesh.txt");
	mesh->Print(geomeshout);
	
	
	TPZCompMesh *cmesh = new TPZCompMesh(mesh);
	
	cmesh->SetDefaultOrder(1);
	cmesh->SetDimModel(3);
	InsertElasticity(cmesh);
	cmesh->AutoBuild();
	
	// Print computatioinal mesh info
	ofstream cmeshout("cmesh.txt");
	cmesh->Print(cmeshout);

	
	TPZSkylineStructMatrix skylstruct(cmesh);
	TPZStepSolver<STATE> step;
	step.SetDirect(ECholesky);
	TPZAnalysis an(cmesh);
	an.SetStructuralMatrix(skylstruct);
	an.SetSolver(step);
	an.Run();
	an.Solution().Print("Solucao");
	
//	TPZAutoPointer<TPZCompMesh> cmeshauto(cmesh);
//	TPZMaterial * mat = cmeshauto->FindMaterial(1);
//	int nstate = mat->NStateVariables();
	int dimension = 2, resolution = 1;
	std::string plotfile("placaaf.vtk");
	TPZVec <std::string> scalnames(2), vecnames(1);
	vecnames[0] = "Displacement";
	scalnames[0] = "Mn1";
	scalnames[1] = "Mn2";
	
	
	an.DefineGraphMesh(dimension, scalnames, vecnames, plotfile);
	an.PostProcess(resolution);

//	TPZMaterial * mat = cmesh->FindMaterial(1);
//	int nstate = mat->NStateVariables();
//	int nscal = 0, nvec = 0, dim = 3;
//	if(nstate ==1) 
//	{
//		nscal = 1;
//	}
//	else
//	{
//		nvec = 1;
//	}
//	TPZManVector<std::string> scalnames(nscal),vecnames(nvec);
//	if(nscal == 1)
//	{
//		scalnames[0]="state";            
//	}
//	else
//	{
//		vecnames[0] = "state";
//	}
//	std::string postprocessname("after.vtk");
//	TPZVTKGraphMesh vtkmesh(cmesh,dim,mat,scalnames,vecnames);
//	vtkmesh.SetFileName(postprocessname);
//	vtkmesh.SetResolution(1);
//	int numcases = 1;
//	
//	
//	// Iteracoes de tempo
//	int istep = 0, nsteps = 2;
//	vtkmesh.DrawMesh(numcases);
//	vtkmesh.DrawSolution(istep, 1.);
	
// // Print geometrical mesh info
//	ofstream geomeshoutf("geomeshf.txt");
//	mesh->Print(geomeshoutf);
//	
//	// Print computatioinal mesh info
//	ofstream cmeshoutf("cmeshf.txt");
//	cmesh->Print(cmeshoutf);
	
	delete cmesh;
	delete mesh;
	return 0;
	
}

TPZGeoMesh *GetMesh (int nx,int ny){
	int i,j;
	int64_t id, index;
	
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
	
	//Auxiliar vector to store a element connectivities
	TPZVec <int64_t> connect(4,0);
	
	//Element connectivities
	for(i = 0; i < (nx - 1); i++){
		for(j = 0; j < (ny - 1); j++){
			index = (i)*(ny - 1)+ (j);
			connect[0] = (i)*ny + (j);
			connect[1] = connect[0]+(ny);
			connect[2] = connect[1]+1;
			connect[3] = connect[0]+1;
			gmesh->CreateGeoElement(EQuadrilateral,connect,1,id);
		}
	}
	//Generate neighborhod information
	gmesh->BuildConnectivity();
	
	int64_t el, numelements = gmesh->NElements();
	int  dirbottID = -1, dirtopID = -2;
	TPZManVector <int64_t> TopolPlate(4);
	
	for (el=0; el<numelements; el++)
	{
		int64_t totalnodes = gmesh->ElementVec()[el]->NNodes();
		TPZGeoEl *plate = gmesh->ElementVec()[el];
		for (int i=0; i<4; i++){
			TopolPlate[i] = plate->NodeIndex(i); 
		}
		
		// Colocando as condicoes de contorno
		TPZManVector <TPZGeoNode> Nodefinder(totalnodes);
		TPZManVector <REAL,3> nodecoord(3);
		// na face x = 1
		TPZVec<int64_t> ncoordzbottVec(0); int64_t sizeOfbottVec = 0;
		TPZVec<int64_t> ncoordztopVec(0); int64_t sizeOftopVec = 0;		
		for (int64_t i = 0; i < totalnodes; i++) 
		{
			Nodefinder[i] = gmesh->NodeVec()[TopolPlate[i]];
			Nodefinder[i].GetCoordinates(nodecoord);
			if (nodecoord[2] == 0. & nodecoord[1] == 0.)
			{
				sizeOfbottVec++;
				ncoordzbottVec.Resize(sizeOfbottVec);
				ncoordzbottVec[sizeOfbottVec-1] = TopolPlate[i];
			}
		    if (nodecoord[2] == 0. & nodecoord[1] == 1.)
			{
				sizeOftopVec++;
				ncoordztopVec.Resize(sizeOftopVec);
				ncoordztopVec[sizeOftopVec-1] = TopolPlate[i];
			}			
		}
		if (sizeOfbottVec == 2) {
			int sidesbott = plate->WhichSide(ncoordzbottVec);
			TPZGeoElSide platesidebott(plate, sidesbott);
			TPZGeoElBC(platesidebott,dirbottID);
		}		
		if (sizeOftopVec == 2) {
			int sidestop = plate->WhichSide(ncoordztopVec);
			TPZGeoElSide platesidetop(plate, sidestop);
			TPZGeoElBC(platesidetop,dirtopID);
		}
		
		ncoordzbottVec.Resize(0);
		sizeOfbottVec = 0;
		ncoordztopVec.Resize(0);
		sizeOftopVec = 0;
	}
	ofstream bf("before.vtk");
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, bf);
	return gmesh;
}

void InsertElasticity(TPZCompMesh *mesh)
{
	mesh->SetDimModel(3);
	int nummat = 1, dirichlet = 0, neumann = 1;
	//int mixed = 2;
	int dirbott = -1, dirtop = -2;
	//TPZManVector<REAL> force(3,0.);
	//	force[2] = 5.;
	REAL h = 0.001;
	REAL f = 0;
	REAL E = 1.7472*10000000;
	REAL ni1 = 0.;
	REAL ni2 = 0.;
	REAL G = E/(2*(1+ni1));
//	TPZFMatrix<REAL> naxes;	
	TPZFMatrix<STATE> naxes(3,3,0.0);
	
	naxes(0,0)=1.0;
	naxes(0,1)=0.0;
	naxes(0,2)=0.0;
	
	naxes(1,0)=0.0;
	naxes(1,1)=1.0;
	naxes(1,2)=0.0;

	naxes(2,0)=0.0;
	naxes(2,1)=0.0;
	naxes(2,2)=1.0;
	
	ofstream file("axes.txt");
	naxes.Print(" axes = ",file, EMathematicaInput);
	
//#ifdef log4cxx
//	if (logdata->isDebugEnabled()) {
//		
//		std::stringstream sout;
//		naxes.Print(" axes ",sout, EMathematicaInput );
//		LOGPZ_DEBUG(lodgata,sout.str());
//	}
//#endif	
	
	TPZVec< STATE > xf(6,0.0);
	//xf[2] = 5;
		
	//TPZMatPoisson3d *p = new TPZMatPoisson3d(nummat, 2);//
    
	//TPZElasticity3D *elast = new TPZElasticity3D(nummat, Ela, poisson, force);//
	
	TPZPlaca *placa = new TPZPlaca(nummat, h, f, E, E, ni1, ni2, G, G, G, naxes, xf);
	
	TPZMaterial * elastauto(placa);
	mesh->InsertMaterialObject(elastauto);
	
	TPZFMatrix<STATE> val1bott(6,6,0.),val2bott(6,1,0.);
	TPZBndCond *bcbott = placa->CreateBC(elastauto, dirbott, dirichlet, val1bott, val2bott);
	TPZMaterial * bcbottauto(bcbott);
	mesh->InsertMaterialObject(bcbottauto);

	TPZFMatrix<STATE> val1top(6,6,0.),val2top(6,1,0.);
	val2top(3.,0.) = 2;
	TPZBndCond *bctop = placa->CreateBC(elastauto, dirtop, neumann, val1top, val2top);
	TPZMaterial * bctopauto(bctop);
	mesh->InsertMaterialObject(bctopauto);
	
	
	/*
	 // Dirichlet em 1 -1 -1 yz;
	 val1(0,0) = 0.;
	 val1(1,1) = 1.;
	 val1(2,2) = 1.;
	 TPZBndCond *bc2 = viscoelast->CreateBC(viscoelastauto, dir2, mixed, val1, val2);
	 TPZMaterial * bcauto2(bc2);
	 mesh->InsertMaterialObject(bcauto2);
	 
	 // Dirichlet em 1 1 -1 z;
	 val1(0,0) = 0.;
	 val1(1,1) = 0.;
	 val1(2,2) = 1.;
	 TPZBndCond *bc3 = viscoelast->CreateBC(viscoelastauto, dir3, mixed, val1, val2);
	 TPZMaterial * bcauto3(bc3);
	 mesh->InsertMaterialObject(bcauto3);
	 */
}

