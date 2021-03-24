/**
 * @file
 * @brief Implements the use of the transform side to side as tutorial example of the integral NeoPZ module
 */
#include <pzvec.h>
#include <pzgmesh.h>
#include <pzquad.h>
#include <pztrnsform.h>

#include "pzgeoel.h"
#include "TPZGeoElement.h"

#include "TPZReadGIDGrid.h"
#include "TPZVTKGeoMesh.h"

#include "TPZRefPatternDataBase.h"

#include "pzlog.h"
//IO
#include <iostream>
#include <algorithm>

using namespace std;

// nx = number of nodes in x direction
// ny = number of nodes in y direction
TPZGeoMesh * GetMesh(int nx,int ny);

void RefineTowards(TPZGeoElSide gelside, int numtimes);

void RefineTowardsCenter(TPZGeoMesh &gmesh, int sidedim, int numtimes);

void RefineTowards(TPZGeoMesh &gmesh, TPZVec<REAL> &x, int sidedim, int numtimes);

std::string DirectorioCorriente = PZSOURCEDIR;

int main() {

    gRefDBase.InitializeRefPatterns();

	DirectorioCorriente += "/Projects/CursoPZ/Class_4/";
/*
	cout << "***********************************************************************\n";
	cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
	cout << "PZ - Class 4 -->> Geometric Meshes and h refinment\n";
	cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
	cout << "***********************************************************************\n\n";
	cout << "Number of nodes: x direction and y direction:     ";
	int nx,ny;
	cin >> nx >> ny;
	
	//Creates the geometric mesh
	TPZGeoMesh *mesh = GetMesh(nx,ny);
	//	mesh->Print(cout);
	TPZVec<TPZGeoEl *> subelindex;
	
	//Divide the first element
	TPZGeoEl *el_0 = mesh->ElementVec()[0];
	if (!el_0) {
		cout << "Null element -->> exit (-1)\n";
		return -1;
	}
	el_0->Divide(subelindex);
	//	mesh->Print(cout);
	//Divide the third subelement
	TPZGeoEl* el_1 = subelindex[2];
	el_1->Divide(subelindex);
	//	mesh->Print(cout);
	//Divide the third subelement
	TPZGeoEl* el_2 = subelindex[2];
	el_2->Divide(subelindex);
	
	//Prints the refined mesh
	mesh->Print(cout);
	
	TPZGeoEl *el_3 = subelindex[2];
	int side = 5;
	int order = 2;
	
	TPZIntPoints *integ = el_3->CreateSideIntegrationRule (side, order);
	int npts = integ->NPoints();
	int i;
	for (i=0;i<npts;i++){
		cout << "=======================================================================" << endl << endl;
		cout << "Ponto de integracao    : " << i << endl << endl;
		cout << "=======================================================================" << endl;
		
		TPZVec<REAL> pt_3(3,0.);
		TPZVec<REAL> pt_0(3,0.);
		TPZVec<REAL> pt_0_aux(3,0.);
		TPZVec<REAL> pt_1(3,0.);
		TPZVec<REAL> pt_2(3,0.);
		TPZVec<REAL> pt_aux(3,0.);
		TPZVec<REAL> pt_aux1(3,0.);
		REAL w = 0.;
		
		// Take the integration point i and its weight
		integ->Point(i, pt_3, w);
		
		//largedim,lowerdim
		TPZTransform<> t_0 (1);
		TPZTransform<> t_0_aux (1);		
		TPZTransform<> t_1 (1);
		TPZTransform<> t_2 (1);
		TPZTransform<> t_3 (1);
		
		//From level 3 to level 0
		cout << "From level 3 to 0 -->> direct operation..." << endl;
		//Take the transformation from side 5 of the subelement to the
		//corresponding side of the specified father (the element in
		//the first level in this case).
		t_0 = el_3->BuildTransform2(5,el_0,t_0);
		//Apply the transformation from the integration point in the
		//subelement rib to a point in rib of the father
		t_0.Apply(pt_3,pt_aux);
		//Take the transformation from the rib to the face...
		t_3 = el_3->SideToSideTransform(5,8);
		t_0 = el_0->SideToSideTransform(5,8);
		t_0.Apply(pt_aux,pt_0);
		t_3.Apply(pt_3,pt_0_aux);
		cout << "Point in master subelement (element level 3).........: " << pt_3 << endl;
		cout << "Point in master subelement after SidetoSideTransform : " << pt_0_aux << endl;
		cout << "Point in master element (element level 0)............: " << pt_0 << endl;
		el_3->X(pt_0_aux,pt_aux);
		cout << "Real point coordinate at subelement..................:" << pt_aux << endl;
		el_0->X(pt_0,pt_aux);
		cout << "Real point coordinate at element.....................:" << pt_aux << endl;
		
		
		//Level by level
		cout << endl << "-----------------------------------------------------------------------" << endl;
		cout << "From level 3 to 0 -->> level by level..." << endl << endl;
		cout << "Level 3 to level 2..." << endl;
		t_2 = el_3->BuildTransform2(5,el_2,t_2);
		t_2.Apply(pt_3,pt_aux1);
		t_2 = el_2->SideToSideTransform(5,8);
		t_2.Apply(pt_aux1,pt_2);
		cout << "Point master element (element level 2)...............: " << pt_2 << endl;
		el_2->X(pt_2,pt_aux);
		cout << "Real coordinate......................................:" << pt_aux << endl << endl;
		pt_2 = pt_aux1;
		
		cout << "Level 2 to level 1..." << endl;		
		t_1 = el_2->BuildTransform2(5,el_1,t_1);
		t_1.Apply(pt_2,pt_aux1);
		t_1 = el_1->SideToSideTransform(5,8);
		t_1.Apply(pt_aux1,pt_1);
		cout << "Point master element (element level 1)...............: " << pt_1 << endl;
		el_1->X(pt_1,pt_aux);
		cout << "Real coordinate......................................:" << pt_aux << endl << endl;
		pt_1 = pt_aux1;
		
		cout << "Level 1 to level 0..." << endl;
		t_0_aux = el_1->BuildTransform2(5,el_0,t_0_aux);
		t_0_aux.Apply(pt_1,pt_aux);
		t_0_aux = el_0->SideToSideTransform(5,8);
		t_0.Apply(pt_aux,pt_0_aux);
		cout << "Point master element (element level 0)...............: " << pt_0_aux << endl;
		el_0->X(pt_0_aux,pt_aux);
		cout << "Real coordinate......................................:" << pt_aux << endl << endl;
		cout << endl << "-----------------------------------------------------------------------" << endl;
	}
*/
	// Malla con GID

	// Ejemplo uni-dimensional para la generacion de una malla para un reservatorio 
	TPZReadGIDGrid grid;
	std::string nombre = DirectorioCorriente;
	nombre += "Reservatorio1D.gid/Reservatorio1D.dump";
	TPZGeoMesh *meshgrid1D = grid.GeometricGIDMesh(nombre);
	if(!meshgrid1D->NElements())
		return 1;
	ofstream meshoutvtk1D("MeshGIDToVTK1D.vtk");
	TPZVTKGeoMesh::PrintGMeshVTK(meshgrid1D,meshoutvtk1D,true);
	meshoutvtk1D.close();

	// Primer ejemplo bi-dimensional para la generacion de una malla para un reservatorio sin material en el topo y el fondo
	nombre = DirectorioCorriente;
	nombre += "Reservatorio2D.gid/Reservatorio2D.dump";
	TPZGeoMesh *meshgrid2D1 = grid.GeometricGIDMesh(nombre);
	if(!meshgrid2D1->NElements())
		return 1;
	ofstream meshoutvtk2D1("MeshGIDToVTK2D1.1.vtk");
	TPZVTKGeoMesh::PrintGMeshVTK(meshgrid2D1,meshoutvtk2D1,true);

	RefineTowardsCenter(*meshgrid2D1, 1, 3);
	ofstream meshoutvtk2D1b("MeshGIDToVTK2D1.2.vtk");
	TPZVTKGeoMesh::PrintGMeshVTK(meshgrid2D1,meshoutvtk2D1b,true);
    
    meshoutvtk2D1.close();
    meshoutvtk2D1b.close();

	// Segundo ejemplo bi-dimensional para la generacion de una malla para un reservatorio con las bordas superior e inferior
	nombre = DirectorioCorriente;
	nombre += "Reservatorio2D2.dump";
	TPZGeoMesh *meshgrid2D2 = grid.GeometricGIDMesh(nombre);
	if(!meshgrid2D2->NElements())
		return 1;
	ofstream meshoutvtk2D2("MeshGIDToVTK2D2.vtk");
	TPZVTKGeoMesh::PrintGMeshVTK(meshgrid2D2,meshoutvtk2D2,true);
	meshoutvtk2D2.close();

	// Ejemplo tri-dimensional para la generacion de una malla para un reservatorio con las bordas superior e inferior
	nombre = DirectorioCorriente;
	nombre += "Reservatorio3D.gid/Reservatorio3D.dump";
	TPZGeoMesh *meshgrid3D = grid.GeometricGIDMesh(nombre);
	if(!meshgrid3D->NElements())
		return 1;
	ofstream meshoutvtk3D("MeshGIDToVTK3D.vtk");
	TPZVTKGeoMesh::PrintGMeshVTK(meshgrid3D,meshoutvtk3D,true);
	meshoutvtk3D.close();

	// Ejemplo de aula
	nombre = DirectorioCorriente;
	nombre += "superfplana.dump";
	TPZGeoMesh *mesh = grid.GeometricGIDMesh(nombre);
	if(!mesh->NElements())
		return 1;
	ofstream meshoutvtk("superfplana.vtk");
	TPZVTKGeoMesh::PrintGMeshVTK(mesh,meshoutvtk,true);
	meshoutvtk.close();
	return 0;
}

TPZGeoMesh *GetMesh (int nx,int ny) {
	int64_t i,j;
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
			coord[2] = coord[0];
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
			//Allocates and define the geometric element
			gmesh->CreateGeoElement(EQuadrilateral,connect,1,id);
		}
	}
	//Generate neighborhod information
	gmesh->BuildConnectivity();
	return gmesh;
}

#include "pzgeoelbc.h"
#include "TPZRefPatternTools.h"

void RefineTowards(TPZGeoElSide gelside, int numtimes)
{
    TPZGeoElBC gbc(gelside,10);
    TPZGeoEl *gel = gelside.Element();
    std::set<int> matids;
    matids.insert(10);
    
    for (int iref=0; iref<numtimes; iref++)
    {
        std::set<TPZGeoEl *> elementstodivide;
        elementstodivide.insert(gel);
        int ncorners = gel->NSideNodes(gelside.Side());
        for (int corner=0; corner<ncorners; corner++) {
            int side = gel->SideNodeLocIndex(gelside.Side(), corner);
            TPZGeoElSide start(gel,side);
            TPZGeoElSide neighbour = start.Neighbour();
            while (neighbour != start) {
                elementstodivide.insert(neighbour.Element());
                neighbour = neighbour.Neighbour();
            }
        }
        std::set<TPZGeoEl *>::iterator it;
        for (it=elementstodivide.begin(); it != elementstodivide.end(); it++) {
            TPZRefPatternTools::RefineDirectional(*it, matids, 11+iref);
        }
    }
}

void RefineTowardsCenter(TPZGeoMesh &gmesh, int sidedim, int numtimes)
{
    if (gmesh.NElements()==0) {
        return;
    }
    TPZManVector<REAL,3> xmin(3), xmax(3),xmid(3);
    
    
    gmesh.NodeVec()[0].GetCoordinates(xmin);
    xmax = xmin;
    int nnodes = gmesh.NNodes();
    for (int in=0; in<nnodes; in++) {
        TPZManVector<REAL,3> xco(3);
        gmesh.NodeVec()[in].GetCoordinates(xco);
        for (int i=0; i<3; i++) {
            if (xmin[i] > xco[i]) {
                xmin[i] = xco[i];
            }
            if (xmax[i] < xco[i]) {
                xmax[i] = xco[i];
            }
        }
    }
    for (int i=0; i<3; i++) {
        xmid[i] = (xmin[i]+xmax[i])/2.;
    }
    RefineTowards(gmesh, xmid, sidedim, numtimes);
}

void RefineTowards(TPZGeoMesh &gmesh, TPZVec<REAL> &x, int sidedim, int numtimes)
{
    if (gmesh.NElements()==0) {
        return;
    }
    int dim = gmesh.ElementVec()[0]->Dimension();
    
    for (int el=0; el<gmesh.NElements(); el++) {
        TPZGeoEl *gel = gmesh.ElementVec()[el];
        if (!gel) {
            continue;
        }
        dim = max(dim,gel->Dimension());
    }

    TPZManVector<REAL,3> qsi(dim,0.);
    int64_t elindex = 0;
    TPZGeoEl *gel = gmesh.FindElement(x, qsi, elindex, dim);
    
    sidedim = min(sidedim,gel->Dimension());
    int side = 0;
    while (gel->SideDimension(side) != sidedim) {
        side++;
    }
    
    TPZGeoElSide targetside(gel,side);
    RefineTowards(targetside, numtimes);
}

