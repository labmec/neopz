#ifdef HAVE_CONFIG_H
#include <pz_config.h>
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

#include "pzskylstrmatrix.h"
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

const int bc1 = -1;
const int bc2 = -2;
const int bc3 = -3;
const int bc4 = -4;

TPZGeoMesh *MalhaGeom(int NRefUnif, REAL Lx, REAL Ly);
TPZCompMesh *MalhaComp(TPZGeoMesh * gmesh,int pOrder);
void BuildHybridMesh(TPZCompMesh *cmesh, std::set<int> &MaterialIDs, int LagrangeMat, int InterfaceMat);
void Solve ( TPZAnalysis &an );



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
//	ofstream arg1("gmesh1.txt");
//	gmesh->Print(arg1);
    	
	TPZCompMesh * cmesh= MalhaComp(gmesh,  p);
	ofstream arg2("cmesh.txt");
	cmesh->Print(arg2);
	
    ofstream arg3("gmesh.txt");
	gmesh->Print(arg3);
    
//    ofstream file("malhageometrica.vtk");
//    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file, true);
    
//    TPZAnalysis an(cmesh);
    
//	Solve( an );


	return EXIT_SUCCESS;
}


TPZGeoMesh *MalhaGeom(int NRefUnif, REAL Lx, REAL Ly)
{
    int64_t Qnodes = 4;
	
	TPZGeoMesh * gmesh = new TPZGeoMesh;
	gmesh->SetMaxNodeId(Qnodes-1);
	gmesh->NodeVec().Resize(Qnodes);
	TPZVec<TPZGeoNode> Node(Qnodes);
	
	TPZVec <int64_t> TopolQuad(4);
	TPZVec <int64_t> TopolLine(2);
	
	//indice dos nos
	int64_t id = 0;
	REAL valx;
	for(int64_t xi = 0; xi < Qnodes/2; xi++)
	{
		valx = xi*Lx;
		Node[id].SetNodeId(id);
		Node[id].SetCoord(0 ,valx );//coord X
		Node[id].SetCoord(1 ,0. );//coord Y
		gmesh->NodeVec()[id] = Node[id];
		id++;
	}
	
	for(int64_t xi = 0; xi < Qnodes/2; xi++)
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
		int64_t n = gmesh->NElements();
		for ( int64_t i = 0; i < n; i++ ){
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
    TPZElasticityHybridMaterial *material = new TPZElasticityHybridMaterial(matInterno, 1., 1., 10., 10.);
    TPZElasticityHybridMaterial *matlagrange = new TPZElasticityHybridMaterial(lagrangemat, 1., 1., 0., 0.);
    TPZElasticityHybridMaterial *matinterface = new TPZElasticityHybridMaterial(interfacemat, 1., 1., 0., 0.);


    
	TPZMaterial * mat1(material);
	TPZMaterial * mat2(matlagrange);
    TPZMaterial * mat3(matinterface);
    
	material->NStateVariables();
	matlagrange->NStateVariables();
   	matinterface->NStateVariables();
    
	TPZCompEl::SetgOrder(pOrder);
	TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDimModel(dim);
    
	cmesh->InsertMaterialObject(mat1);
	cmesh->InsertMaterialObject(mat2);
    cmesh->InsertMaterialObject(mat3);
	
	///Inserir condicao de contorno
	TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
	REAL uN=0.;
	val2(0,0)=uN;
	TPZMaterial * BCondN1 = material->CreateBC(mat1, bc1,neumann, val1, val2);
	cmesh->InsertMaterialObject(BCondN1);
    
    TPZMaterial * BCondN2 = material->CreateBC(mat1, bc3,neumann, val1, val2);
    cmesh->InsertMaterialObject(BCondN2);
	
	TPZFMatrix<STATE> val12(2,2,0.), val22(2,1,0.);
	REAL uD=0.;
	val22(0,0)=uD;
	TPZMaterial * BCondD1 = material->CreateBC(mat1, bc2,dirichlet, val12, val22);
	cmesh->InsertMaterialObject(BCondD1);
	
	TPZMaterial * BCondD2 = material->CreateBC(mat1, bc4,dirichlet, val12, val22);
	cmesh->InsertMaterialObject(BCondD2);
	
    cmesh->SetAllCreateFunctionsContinuous();

	set<int> SETmat1;
	SETmat1.insert(bc1);
    SETmat1.insert(bc2);
    SETmat1.insert(bc3);
    SETmat1.insert(bc4);
    
	//criar set dos materiais
    std::set<int> MaterialIDs;
    std::set<int> BCMaterialIDs;
    MaterialIDs.insert(matInterno);
//     MaterialIDs.insert(lagrangemat);
//     MaterialIDs.insert(interfacemat);
//     MaterialIDs.insert(bc1);
//     MaterialIDs.insert(bc2);
//     MaterialIDs.insert(bc3);
//     MaterialIDs.insert(bc4);
    BCMaterialIDs.insert(bc1);
    BCMaterialIDs.insert(bc2);
    BCMaterialIDs.insert(bc3);
    BCMaterialIDs.insert(bc4);

    
    TPZBuildMultiphysicsMesh::BuildHybridMesh(cmesh, MaterialIDs, BCMaterialIDs, lagrangemat, interfacemat);
        
	return cmesh;
}

void Solve ( TPZAnalysis &an ){
	TPZCompMesh *malha = an.Mesh();
    
    //TPZBandStructMatrix mat(malha);
	//TPZSkylineStructMatrix mat(malha);// requer decomposição simétrica, não pode ser LU!
	//TPZBlockDiagonalStructMatrix mat(malha);//ok!
	//TPZFrontStructMatrix<TPZFrontNonSym> mat ( malha );// não funciona com método iterativo
	TPZFStructMatrix mat( malha );// ok! matriz estrutural cheia
	//TPZSpStructMatrix mat( malha );//matriz estrutural esparsa (???? NÃO FUNCIONOU !!!!!!!!!!)
	TPZStepSolver<STATE> solv;
	solv.SetDirect (  ELU );//ECholesky);// ELU , ELDLt ,
    
	
    //	cout << "ELDLt " << endl;
	an.SetSolver ( solv );
	an.SetStructuralMatrix ( mat );
	cout << endl;
	an.Solution().Redim ( 0,0 );
	cout << "Assemble " << endl;
	an.Assemble();
    //	std::ofstream fileout("rigidez.txt");
    //	an.Solver().Matrix()->Print("Rigidez", fileout, EMathematicaInput);
	
//	an.Solve();
//	cout << endl;
//	cout << "No equacoes = " << malha->NEquations() << endl;
}

