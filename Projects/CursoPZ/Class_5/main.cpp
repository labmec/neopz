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
#include "pzbndcond.h"

#include "pzfstrmatrix.h"
#include "TPZPersistenceManager.h"

#include <iostream>
#include <fstream>
using namespace std;

// nx = number of nodes in x direction
// ny = number of nodes in y direction
TPZGeoMesh * GetMesh(int nx,int ny);

int main() {
	//TPZSavable::Register(TPZSAVEABLEID,Restore<TPZSavable>);
	cout << "***********************************************************************\n";
	cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
	cout << "PZ - Class 5 -->> Writing and reading meshes on files\n";
	cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
	cout << "***********************************************************************\n\n";
	cout << "Number of nodes: x direction and y direction:     ";
	int nx,ny;
	//cin >> nx >> ny;
    nx =3;
    ny = 3;
	int order;
	cout << "Enter the order of the element:\n";
//	cin >> order;
    order = 1;
	TPZCompEl::SetgOrder(order);
	
	//Creates the geometric mesh
	TPZGeoMesh *gmesh = GetMesh(nx,ny);
    gmesh->ElementVec()[0]->CreateBCGeoEl(4, -1);
    gmesh->ElementVec()[2]->CreateBCGeoEl(4, -2);
    gmesh->ElementVec()[2]->CreateBCGeoEl(5, -3);
    gmesh->ElementVec()[3]->CreateBCGeoEl(2, -4);
	gmesh->SetName("testing a space");
	
	TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(2);
	
	TPZMat2dLin *mat2d = new TPZMat2dLin (1);
	TPZFMatrix<STATE> xkin (1,1,1);
	TPZFMatrix<STATE> xcin (1,1,1.);
	TPZFMatrix<STATE> xfin (1,1,1e3);
	mat2d->SetMaterial(xkin,xcin,xfin);
	cmesh->InsertMaterialObject (mat2d);
    
    TPZFMatrix<STATE> val1(1,1,0.),val2(1,1,1.);
    TPZBndCond *bc1 = mat2d->CreateBC(mat2d, -1, 1, val1, val2);
    cmesh->InsertMaterialObject(bc1);
    val1(0,0) = 1.;
    val2(0,0) = 5.;
    TPZBndCond *bc2 = mat2d->CreateBC(mat2d, -2, 2, val1, val2);
    cmesh->InsertMaterialObject(bc2);
    
    val2(0,0) = 5;
    TPZBndCond *bc3 = mat2d->CreateBC(mat2d, -3, 0, val1, val2);
    cmesh->InsertMaterialObject(bc3);
    
    val2(0,0) = 0.;
    TPZBndCond *bc4 = mat2d->CreateBC(mat2d, -4, 0, val1, val2);
    cmesh->InsertMaterialObject(bc4);
	//TPZMaterial* mat(mat2d);
//    cmesh->SetAllCreateFunctionsDiscontinuous();
//    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->SetAllCreateFunctionsHDivPressure();
	cmesh->AutoBuild();
//    cmesh->ApproxSpace().CreateInterfaces(*cmesh);
	TPZVec<int64_t> subelindex(4,0);
    int64_t elindex = 0;
    int shouldinterpolate = 0;
	cmesh->ElementVec()[elindex]->Divide(elindex,subelindex,shouldinterpolate);

    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    int64_t neq = cmesh->NEquations();
    TPZFMatrix<STATE> rhs(neq,1);
    TPZFStructMatrix full(cmesh);
    TPZMatrix<STATE> *stiff = full.CreateAssemble(rhs, 0);

    
    {
        ofstream out("all.dat");
        cmesh->Reference()->Print(out);
        out << "antes de lido\n";
        cmesh->Print(out);
        stiff->Print("Rigidez global",out);
    }
	{
                TPZPersistenceManager::OpenWrite("dump.dat");
                TPZPersistenceManager::WriteToFile(gmesh);
                TPZPersistenceManager::WriteToFile(cmesh);
                TPZPersistenceManager::CloseWrite();
                
	}
	{
                TPZPersistenceManager::OpenRead("dump.dat");
		TPZGeoMesh *tst = dynamic_cast<TPZGeoMesh*>(TPZPersistenceManager::ReadFromFile());
		TPZCompMesh *tsc = dynamic_cast<TPZCompMesh *>(TPZPersistenceManager::ReadFromFile());
        std::cout << "depois de lido "<<endl;
		//if(tst) tst->Print(out);
		if(tsc) tsc->Print(std::cout);
		delete tst;
	}
	if(cmesh) delete cmesh;
	if(gmesh) delete gmesh;
	return 0;
	
}

TPZGeoMesh *GetMesh (int nx,int ny){
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
	return gmesh;
}
