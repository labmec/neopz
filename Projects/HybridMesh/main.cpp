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
#include "tpzautopointer.h"
#include "pzbndcond.h"
#include "pzanalysis.h"

#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include "TPZStructMatrix.h"
#include "pzfstrmatrix.h"
#include "TPZFrontNonSym.h"
#include "TPZFrontSym.h"
#include "TPZBSpStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzbstrmatrix.h"

#include "pzpoisson3d.h"
#include "pzhybridpoisson.h"
#include "pzpoisson3dreferred.h"
#include "pzelasmat.h" 
#include "pzelasthybrid.h"

#include "pzbuildmultiphysicsmesh.h"
#include "pzelementgroup.h"
#include "pzcondensedcompel.h"

#include "pzlog.h"

#include "TPZVTKGeoMesh.h"

#include <iostream>
#include <string>

#include <math.h>
#include <set>

#ifdef PZ_LOG
static TPZLogger logger("pz.multiphysics");
#endif

using namespace std;

const int matInterno = 1;
const int lagrangemat = 2;
const int interfacemat = 3;

const int dirichlet = 0;
//const int neumann = 1;
//const int mixed = 2;

const int bc1 = -1;
const int bc2 = -2;
const int bc3 = -3;
const int bc4 = -4;

TPZGeoMesh *MalhaGeom(int NRefUnif, REAL Lx, REAL Ly);
TPZCompMesh *MalhaComp(TPZGeoMesh * gmesh,int pOrder);

void BuildElementGroups(TPZCompMesh *cmesh, int materialid, int interfacemat, int lagrangemat);

void ResetMesh(TPZCompMesh *cmesh);

int main(int argc, char *argv[])
{
	
	int  p=1;
	int  NRefUnif=1;
    REAL Lx=1.;
    REAL Ly=1.;
	
    TPZGeoMesh * gmesh = MalhaGeom(NRefUnif,Lx,Ly);
	ofstream arg1("gmesh1.txt");
	gmesh->Print(arg1);
    	
	TPZCompMesh * cmesh= MalhaComp(gmesh, p);

    int neq = cmesh->NEquations();
    TPZFMatrix<STATE> stiff(neq,neq,0.),rhs(neq,1,0.);
    TPZFStructMatrix<STATE> fstr(cmesh);
    fstr.Assemble(stiff, rhs, 0);
    
    ofstream arg4("stiffness.txt");
    stiff.Print("Global Stiffness matrix",arg4);
	
    TPZAnalysis an(cmesh);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    TPZFStructMatrix<STATE> fullstr(cmesh);
    an.SetStructuralMatrix(fullstr);
    an.SetSolver(step);
    an.Run();

    ofstream arg2("cmesh.txt");
	cmesh->Print(arg2);

    
    ofstream arg3("gmesh.txt");
	gmesh->Print(arg3);
    
    ofstream file("malhageometrica.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file, true);
    
    ResetMesh(cmesh);
    delete cmesh;
    delete gmesh;
    
	return EXIT_SUCCESS;
}

TPZGeoMesh *MalhaGeom(int NRefUnif, REAL Lx, REAL Ly)
{
    int Qnodes = 4;
	
	TPZGeoMesh * gmesh = new TPZGeoMesh;
	gmesh->SetMaxNodeId(Qnodes-1);
	gmesh->NodeVec().Resize(Qnodes);
	TPZVec<TPZGeoNode> Node(Qnodes);
	
	TPZVec <int64_t> TopolQuad(4);
	TPZVec <int64_t> TopolLine(2);
	
	//indice dos nos
	int64_t id = 0;
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
    gmesh->AddInterfaceMaterial(matInterno,lagrangemat, interfacemat);
//    gmesh->AddInterfaceMaterial(lagrangemat,bc1, bc1);
//    gmesh->AddInterfaceMaterial(lagrangemat,bc2, bc2);
//    gmesh->AddInterfaceMaterial(lagrangemat,bc3, bc3);
//    gmesh->AddInterfaceMaterial(lagrangemat,bc4, bc4);

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

TPZCompMesh* MalhaComp(TPZGeoMesh * gmesh, int pOrder)
{
	/// criar materiais
	int dim = 2;
    TPZMatPoisson3d *material = new TPZMatPoisson3d(matInterno,dim);
    TPZMatPoisson3d *matlagrange = new TPZHybridPoisson(lagrangemat, dim);
    TPZMatPoisson3d *matinterface = new TPZHybridPoisson(interfacemat,dim);
    
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
	TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
	
	TPZFMatrix<STATE> val12(2,2,0.), val22(2,1,0.);
	REAL uD=0.;
	val22(0,0)=uD;
	TPZMaterial * BCondD1 = material->CreateBC(mat1, bc2,dirichlet, val12, val22);
	cmesh->InsertMaterialObject(BCondD1);
	
	TPZMaterial * BCondD2 = material->CreateBC(mat1, bc4,dirichlet, val12, val22);
	cmesh->InsertMaterialObject(BCondD2);
	
	REAL uN=0.;
	val2(0,0)=uN;
    
	TPZMaterial * BCondN1 = material->CreateBC(mat1, bc1,dirichlet, val1, val2);
	cmesh->InsertMaterialObject(BCondN1);
    
    TPZMaterial * BCondN2 = material->CreateBC(mat1, bc3,dirichlet, val1, val2);
    cmesh->InsertMaterialObject(BCondN2);
    
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
    BuildElementGroups(cmesh, matInterno, interfacemat,lagrangemat);
    
    //    TPZMaterial *lagrange = cmesh->MaterialVec()[lagrangemat];
    //    cmesh->MaterialVec().erase(lagrangemat);
    //    delete lagrange;
    
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

void BuildElementGroups(TPZCompMesh *cmesh, int materialid, int interfacemat, int lagrangemat)
{
    int64_t nel = cmesh->NElements();
    std::map<int,TPZElementGroup *> elgroup;
    cmesh->LoadReferences();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = cmesh->ElementVec()[el];
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (gel && gel->MaterialId() == materialid) {
            int64_t index;
            TPZElementGroup *elgr = new TPZElementGroup(*cmesh,index);
            elgroup[el] = elgr;
#ifdef PZ_LOG
            {
                std::stringstream sout;
                sout << "Creating an element group around element index " << el;
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif
        }
    }
    for (std::map<int,TPZElementGroup *>::iterator it = elgroup.begin(); it != elgroup.end(); it++) {
        TPZCompEl *cel = cmesh->ElementVec()[it->first];
        TPZGeoEl *gel = cel->Reference();
        int dim = gel->Dimension();
        int nsides = gel->NSides();
        for (int is=0; is<nsides; is++) {
            int sidedim = gel->SideDimension(is);
            if (sidedim != dim-1) {
                continue;
            }
            TPZStack<TPZCompElSide> equal;
            TPZGeoElSide gelside(gel,is);
            gelside.EqualLevelCompElementList(equal, 0, 0);
            // look for the interface elements (which look at me)
            int64_t neq = equal.size();
            for (int64_t eq =0; eq < neq; eq++) {
                TPZCompElSide eqsize = equal[eq];
                TPZCompEl *eqcel = eqsize.Element();
                TPZGeoEl *eqgel = eqcel->Reference();
                if (!eqgel) {
                    continue;
                }
                int matid = eqgel->MaterialId();
                if (matid == interfacemat) {
                    TPZInterfaceElement *interfacee = dynamic_cast<TPZInterfaceElement *>(eqcel);
                    if (!interfacee ) {
                        DebugStop();
                    }
                    TPZCompEl *left = interfacee->LeftElement();
                    TPZCompEl *right = interfacee->RightElement();
                    if (left == cel || right == cel) {
                        it->second->AddElement(interfacee);
                        eqgel->ResetReference();
                    }
                }
                // look for the boundary elements (all elements which are not lagrange materials)
                else if (matid != lagrangemat && eqgel->Dimension() == dim-1)
                {
                    it->second->AddElement(eqcel);
                    eqgel->ResetReference();
                }
            }
        }
        it->second->AddElement(cel);
        gel->ResetReference();
    }
    cmesh->ComputeNodElCon();
    nel = cmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = cmesh->ElementVec()[el];
        TPZElementGroup *elgr = dynamic_cast<TPZElementGroup *>(cel);
        if(elgr) {
            TPZCondensedCompEl *cond = NULL;
            cond = new TPZCondensedCompEl(elgr);
        }
    }
    cmesh->CleanUpUnconnectedNodes();
}

void ResetMesh(TPZCompMesh *cmesh)
{
    int64_t nel = cmesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = cmesh->ElementVec()[el];
        TPZCondensedCompEl *cond = dynamic_cast<TPZCondensedCompEl *>(cel);
        if (cond) {
            cond->Unwrap();
        }
    }
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = cmesh->ElementVec()[el];
        TPZElementGroup *elgr = dynamic_cast<TPZElementGroup *>(cel);
        if (elgr) {
            elgr->Unwrap();
        }
    }
}
