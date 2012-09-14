#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "pzvec.h"
#include "pzstack.h"
#include "pzfmatrix.h"
#include "pzfstrmatrix.h"
#include "pzlog.h"

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "TPZInterfaceEl.h"
#include "pzgeoelside.h"
#include "TPZGeoLinear.h"
#include "pzgeopoint.h"

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "tpzcompmeshreferred.h"
#include "tpzautopointer.h"
#include "pzbndcond.h"
#include "pzanalysis.h"

#include "TPZParSkylineStructMatrix.h"
#include "pzstepsolver.h"
#include "pzstrmatrix.h"
#include "TPZFrontNonSym.h"
#include "TPZFrontSym.h"
#include "TPBSpStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzbstrmatrix.h"

#include "pzpoisson3d.h"
#include "pzpoisson3dreferred.h"

#include "pzbuildmultiphysicsmesh.h"

#include "pzlog.h"

#include <iostream>
#include <string>

#include <math.h>
#include <set>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.multiphysics"));
#endif

using namespace std;

const int matInterno = 1;
const int lagrangemat = 2;
const int interfacemat = 3;

const int dirichlet = 0;
const int neumann = 1;
const int mixed = 2;

const int bc1 = -1;
const int bc2 = -2;
const int bc3 = -3;
const int bc4 = -4;

TPZGeoMesh *MalhaGeom(bool interface1, int nh, REAL w, REAL L);
TPZCompMesh *MalhaCompComInterf(TPZGeoMesh * gmesh,int pOrder);
void BuildHybridMesh(TPZCompMesh *cmesh, std::set<int> &MaterialIDs, int LagrangeMat, int InterfaceMat);

void PrintGMeshVTK(TPZGeoMesh * gmesh, std::ofstream &file);


int main(int argc, char *argv[])
{
#ifdef LOG4CXX
	std::string logs("log4cxx.stabilizedhybrid");
	InitializePZLOG();
#endif
	
	int p =1;
	
	// -------- validar o metodo BuildHybridMesh() ------------
	TPZGeoMesh * gmesh = MalhaGeom(true,1,1.,1.);
	ofstream arg1("mygmesh.txt");
	gmesh->Print(arg1);
    
    ofstream file1("malhageoRefUniforme.vtk");
	PrintGMeshVTK(gmesh, file1);
	
	TPZCompMesh * cmesh= MalhaCompComInterf(gmesh,  p);
	ofstream arg2("mycmesh.txt");
	cmesh->Print(arg2);
	
    ofstream arg12("mygmesh2.txt");
	gmesh->Print(arg12);

	return EXIT_SUCCESS;
}


TPZGeoMesh *MalhaGeom(bool interface1, int nh, REAL w, REAL L)
{
    int Qnodes = 4;
	
	TPZGeoMesh * gmesh = new TPZGeoMesh;
	gmesh->SetMaxNodeId(Qnodes-1);
	gmesh->NodeVec().Resize(Qnodes);
	TPZVec<TPZGeoNode> Node(Qnodes);
	
	TPZVec <int> TopolQuad(4);
	TPZVec <int> TopolLine(2);
	
	//indice dos nos
	int id = 0;
	REAL valx;
	for(int xi = 0; xi < Qnodes/2; xi++)
	{
		valx = xi*L;
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0 ,valx );//coord X
		Node[id].SetCoord(1 ,0. );//coord Y
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}
	
	for(int xi = 0; xi < Qnodes/2; xi++)
	{
		valx = L - xi*L;
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0 ,valx );//coord X
		Node[id].SetCoord(1 ,w);//coord Y
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}

	//indice dos elementos
	id = 0;
    
    //elementos internos
    TopolQuad[0] = 0;
	TopolQuad[1] = 1;
	TopolQuad[2] = 2;
	TopolQuad[3] = 3;
	new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (id,TopolQuad,matInterno,*gmesh);
	id++;
	    
    //elementos de contorno
	TopolLine[0] = 0;
	TopolLine[1] = 1;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc1,*gmesh);
	id++;
	
	TopolLine[0] = 1;
	TopolLine[1] = 2;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc2,*gmesh);
	id++;
	
	TopolLine[0] = 2;
	TopolLine[1] = 3;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc3,*gmesh);
	id++;
	
	TopolLine[0] = 3;
	TopolLine[1] = 0;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,bc4,*gmesh);
	id++;
    
    //construir a malha
	gmesh->BuildConnectivity();
	
    
    //Refinamento uniforme
	for( int ref = 0; ref < nh; ref++ ){
		TPZVec<TPZGeoEl *> filhos;
		int n = gmesh->NElements();
		for ( int i = 0; i < n; i++ ){
			TPZGeoEl * gel = gmesh->ElementVec()[i];
            gel->Divide (filhos);
		}//for i
	}//ref
		
	return gmesh;
	
}


TPZCompMesh*MalhaCompComInterf(TPZGeoMesh * gmesh, int pOrder)
{
	/// criar materiais
	int dim = 2;
	
	TPZMatPoisson3d *material = new TPZMatPoisson3d(matInterno,dim);
	TPZMatPoisson3d *matlagrange = new TPZMatPoisson3d(lagrangemat,1); 
    TPZMatPoisson3d *matinterface = new TPZMatPoisson3d(interfacemat,1); 
	
	TPZMaterial * mat1(material);
	TPZMaterial * mat2(matlagrange);
    TPZMaterial * mat3(matinterface);
    
	REAL diff = -1.;
	REAL conv = 0.;
	TPZVec<REAL> convdir(3,0.);
	REAL flux = 8.;
	
	material->SetParameters(diff, conv, convdir);
	material->SetInternalFlux( flux);
	material->NStateVariables();
	
	TPZCompEl::SetgOrder(pOrder);
	TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDimModel(dim);
    
	cmesh->InsertMaterialObject(mat1);
	cmesh->InsertMaterialObject(mat2);
    cmesh->InsertMaterialObject(mat3);
	
	///Inserir condicao de contorno
	TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
	REAL uN=0.;
	val2(0,0)=uN;
	TPZMaterial * BCondN1 = material->CreateBC(mat1, bc1,neumann, val1, val2);
	cmesh->InsertMaterialObject(BCondN1);
    
    TPZMaterial * BCondN2 = material->CreateBC(mat1, bc3,neumann, val1, val2);
    cmesh->InsertMaterialObject(BCondN2);
	
	TPZFMatrix<REAL> val12(2,2,0.), val22(2,1,0.);
	REAL uD=0.;
	val22(0,0)=uD;
	TPZMaterial * BCondD1 = material->CreateBC(mat1, bc2,dirichlet, val12, val22);
	cmesh->InsertMaterialObject(BCondD1);
	
	TPZMaterial * BCondD2 = material->CreateBC(mat1, bc4,dirichlet, val12, val22);
	cmesh->InsertMaterialObject(BCondD2);
	
    cmesh->SetAllCreateFunctionsContinuous();

    //AQUI: AutoBuild(matId)
	set<int> SETmat1;
	SETmat1.insert(bc1);
    SETmat1.insert(bc2);
    SETmat1.insert(bc3);
    SETmat1.insert(bc4);
	cmesh->AutoBuild(SETmat1);
	gmesh->ResetReference();
    
	//criar elementos de interface e 1d (lagrange)
    std::set<int> MaterialIDs;
    MaterialIDs.insert(matInterno);
    MaterialIDs.insert(lagrangemat);
    MaterialIDs.insert(interfacemat);
    
//    set<int>::iterator it;
//    cout << "myset contains:";
//    for (it=MaterialIDs.begin() ; it != MaterialIDs.end(); it++)
//        cout << " " << *it;
//        cout << endl;
    
    TPZBuildMultiphysicsMesh::BuildHybridMesh(cmesh, MaterialIDs, lagrangemat, interfacemat);
    
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
	return cmesh;
}

void PrintGMeshVTK(TPZGeoMesh * gmesh, std::ofstream &file)
{
	file.clear();
	int nelements = gmesh->NElements();
	
	std::stringstream node, connectivity, type;
	
	//Header
	file << "# vtk DataFile Version 3.0" << std::endl;
	file << "TPZGeoMesh VTK Visualization" << std::endl;
	file << "ASCII" << std::endl << std::endl;
	
	file << "DATASET UNSTRUCTURED_GRID" << std::endl;
	file << "POINTS ";
	
	int actualNode = -1, size = 0, nVALIDelements = 0;
	
	for(int el = 0; el < nelements; el++)
	{        
		if(gmesh->ElementVec()[el]->Type() == EPoint)//Exclude Lines and Arc3D
		{
			continue;
		}
		if(gmesh->ElementVec()[el]->Type() == EOned)//Exclude Lines and Arc3D
		{
			continue;
		}
		if(gmesh->ElementVec()[el]->HasSubElement())
		{
			continue;
		}
		
		int elNnodes = gmesh->ElementVec()[el]->NNodes();
		size += (1+elNnodes);
		connectivity << elNnodes;
		
		for(int t = 0; t < elNnodes; t++)
		{
			for(int c = 0; c < 3; c++)
			{
				double coord = gmesh->NodeVec()[gmesh->ElementVec()[el]->NodeIndex(t)].Coord(c);
				node << coord << " ";
			}            
			node << std::endl;
			
			actualNode++;
			connectivity << " " << actualNode;
		}
		connectivity << std::endl;
		
		int elType = -1;
		switch (gmesh->ElementVec()[el]->Type())
		{
			case (ETriangle):
			{
				elType = 5;
				break;                
			}
			case (EQuadrilateral ):
			{
				elType = 9;
				break;                
			}
			case (ETetraedro):
			{
				elType = 10;
				break;                
			}
			case (EPiramide):
			{
				elType = 14;
				break;                
			}
			case (EPrisma):
			{
				elType = 13;
				break;                
			}
			case (ECube):
			{
				elType = 12;
				break;                
			}
			default:
			{
				//ElementType NOT Found!!!
				DebugStop();
				break;    
			}
		}
		
		type << elType << std::endl;
		nVALIDelements++;
	}
	node << std::endl;
	actualNode++;
	file << actualNode << " float" << std::endl << node.str();
	
	file << "CELLS " << nVALIDelements << " ";
	
	file << size << std::endl;
	file << connectivity.str() << std::endl;
	
	file << "CELL_TYPES " << nVALIDelements << std::endl;
	file << type.str();
	
	file.close();
}