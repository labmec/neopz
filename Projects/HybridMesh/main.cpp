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
#include "pzelasmat.h" 
#include "pzelasthybrid.h"

#include "pzbuildmultiphysicsmesh.h"

#include "pzlog.h"

#include "TPZVTKGeoMesh.h"

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

TPZGeoMesh *MalhaGeom(int NRefUnif, REAL Lx, REAL Ly);
TPZCompMesh *MalhaComp(TPZGeoMesh * gmesh,int pOrder);


int main(int argc, char *argv[])
{
#ifdef LOG4CXX
	InitializePZLOG();
#endif
	
	int  p=1;
	int  NRefUnif=1;
    REAL Lx=1.;
    REAL Ly=1.;
	
    TPZGeoMesh * gmesh = MalhaGeom(NRefUnif,Lx,Ly);
	ofstream arg1("gmesh1.txt");
	gmesh->Print(arg1);
    	
	TPZCompMesh * cmesh= MalhaComp(gmesh, p);
	ofstream arg2("cmesh.txt");
	cmesh->Print(arg2);
	
    ofstream arg3("gmesh.txt");
	gmesh->Print(arg3);
    
    ofstream file("malhageometrica.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file, true);
    
	return EXIT_SUCCESS;
}


TPZGeoMesh *MalhaGeom(int NRefUnif, REAL Lx, REAL Ly)
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
		valx = xi*Lx;
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0 ,valx );//coord X
		Node[id].SetCoord(1 ,0. );//coord Y
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}
	
	for(int xi = 0; xi < Qnodes/2; xi++)
	{
		valx = Lx - xi*Lx;
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0 ,valx );//coord X
		Node[id].SetCoord(1 ,Ly);//coord Y
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
    
    // define entre quais materiais vou criar interfaces e o terceiro argumento é o tipo de material que quero nessa interface.
    gmesh->AddInterfaceMaterial(lagrangemat,matInterno, interfacemat);
    gmesh->AddInterfaceMaterial(lagrangemat,bc1, bc1);
    gmesh->AddInterfaceMaterial(lagrangemat,bc2, bc2);
    gmesh->AddInterfaceMaterial(lagrangemat,bc3, bc3);
    gmesh->AddInterfaceMaterial(lagrangemat,bc4, bc4);

    //construir a malha
	gmesh->BuildConnectivity();
	
    
    //Refinamento uniforme
	for( int ref = 0; ref < NRefUnif; ref++ ){
		TPZVec<TPZGeoEl *> filhos;
		int n = gmesh->NElements();
		for ( int i = 0; i < n; i++ ){
			TPZGeoEl * gel = gmesh->ElementVec()[i];
            gel->Divide (filhos);
		}//for i
	}//ref
		
	return gmesh;
	
}


TPZCompMesh*MalhaComp(TPZGeoMesh * gmesh, int pOrder)
{
	/// criar materiais
	int dim = 2;
    TPZMatPoisson3d *material = new TPZMatPoisson3d(matInterno,dim);
    TPZMatPoisson3d *matlagrange = new TPZMatPoisson3d(lagrangemat, dim);
    TPZMatPoisson3d *matinterface = new TPZMatPoisson3d(interfacemat,dim);
    
	TPZMaterial * mat1(material);
	TPZMaterial * mat2(matlagrange);
    TPZMaterial * mat3(matinterface);
    
	material->NStateVariables();
	matlagrange->NStateVariables();
   	matinterface->NStateVariables();
    
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
	
    cmesh->SetAllCreateFunctionsHDivPressure();

	set<int> SETmat1;
	SETmat1.insert(bc1);
    SETmat1.insert(bc2);
    SETmat1.insert(bc3);
    SETmat1.insert(bc4);
    
	//criar set dos materiais
    std::set<int> MaterialIDs;
    std::set<int> BCMaterialIDs;
    MaterialIDs.insert(matInterno);
    //MaterialIDs.insert(lagrangemat);
    //MaterialIDs.insert(interfacemat);
    BCMaterialIDs.insert(bc1);
    BCMaterialIDs.insert(bc2);
    BCMaterialIDs.insert(bc3);
    BCMaterialIDs.insert(bc4);

//    cmesh->AutoBuild(MaterialIDs);
    
    TPZBuildMultiphysicsMesh::BuildHybridMesh(cmesh, MaterialIDs, BCMaterialIDs, lagrangemat, interfacemat);
    
    
//    int nel = cmesh->NElements();
//    for(int i=0; i<nel; i++){
//        TPZCompEl *cel = cmesh->ElementVec()[i];
//        if(!cel) continue;
//        
//        int mid = cel->Material()->Id();
//       
//        if(mid==lagrangemat){
//            
//            int nsides = cel->Reference()->NSides();
//            
//            for(int i = 0; i<nsides; i++){
//                TPZConnect &newnod = cel->Connect(i);
//                newnod.SetPressure(true);
//            }
//        }
//    }
    
    return cmesh;
}
