/**
 * @file
 * @brief Implements the use of the transform side to side as tutorial example of the integral NeoPZ module
 */
#include <pzvec.h>
#include <pzgmesh.h>
#include <pzgeoel.h>
#include <pzquad.h>
#include <pztrnsform.h>
#include <TPZGeoElement.h>


//IO
#include <iostream>
using namespace std;

// nx = number of nodes in x direction
// ny = number of nodes in y direction
TPZGeoMesh * GetMesh(int nx,int ny);

int main(){
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
	TPZGeoEl * el_1 = subelindex[2];
	el_1->Divide(subelindex);
//	mesh->Print(cout);
	//Divide the third subelement
	TPZGeoEl * el_2 = subelindex[2];
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
		TPZTransform t_0 (1);
		TPZTransform t_0_aux (1);		
		TPZTransform t_1 (1);
  	TPZTransform t_2 (1);
	  TPZTransform t_3 (1);

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
	return 0;
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
			coord[2] = coord[0];
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
  return gmesh;
}
