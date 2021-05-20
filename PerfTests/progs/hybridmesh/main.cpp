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
#include "TPZLinearAnalysis.h"

#include "pzstepsolver.h"
#include "pzstrmatrix.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
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


#include "run_stats_table.h"
#include "arglib.h"

RunStatsTable total_rdt("-tot_rdt","Statistics for the whole application");
RunStatsTable cond_rdt("-cond_rdt","Statistics for the static condensation step");
RunStatsTable cond_ass_rdt("-cond_ass_rdt","Statistics for the assemble step during static condensation");

clarg::argBool help("-h",  "display the help message");
clarg::argBool cond_f("-cond",  "perform static condensation");
clarg::argInt  p_order("-p", "polynomial order",1);
clarg::argInt  n_uref("-nuref", "Number of uniform refinements",1);
clarg::argInt  n_threads("-nthreads", "Number of threads",1);

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

REAL const Pi = 4.*atan(1.);

TPZGeoMesh *MalhaGeom(int NRefUnif, REAL Lx, REAL Ly);
TPZCompMesh *MalhaComp(TPZGeoMesh * gmesh,int pOrder);

//funcao f(x)
void Forcing(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);

//sol exata
void SolExata(const TPZVec<REAL> &pt, TPZVec<STATE> &disp,TPZFMatrix<STATE> &flux);

void PosProcessSol(TPZLinearAnalysis &an, std::string plotfile);

void BuildElementGroups(TPZCompMesh *cmesh, int materialid, int interfacemat, int lagrangemat);

void ResetMesh(TPZCompMesh *cmesh);


int main(int argc, char *argv[])
{
    clarg::parse_arguments(argc, argv);
    
    if (help.get_value() == true) {
        cout << "Usage: " << argv[0] << endl;
        clarg::arguments_descriptions(cout, "  ", "\n");
        return 1;
    }
 
	
    REAL Lx=1.;
    REAL Ly=1.;
	
    TPZGeoMesh * gmesh = MalhaGeom(n_uref.get_value(),Lx,Ly);

#ifdef PZDEBUG
    std::cout << "Number of elements = " << gmesh->NElements() << std::endl;
	ofstream arg1("gmesh0.txt");
	gmesh->Print(arg1);
#endif
    	
	TPZCompMesh * cmesh= MalhaComp(gmesh, p_order.get_value());
    
#ifdef PZDEBUG
    std::cout << "Number of equations = " << cmesh->NEquations() << std::endl;
#endif
    
    total_rdt.start();

//    if (cond_f.was_set()) {
        // Perform static condensation
        cond_rdt.start();
        BuildElementGroups(cmesh, matInterno, interfacemat,lagrangemat);
        int neq = cmesh->NEquations();
    
        TPZFMatrix<STATE> stiff(neq,neq,0.),rhs(neq,1,0.);
        TPZSkylineStructMatrix fstr(cmesh);
        cond_ass_rdt.start();
        fstr.SetNumThreads(n_threads.get_value());
        fstr.Assemble(stiff, rhs, 0);
        cond_ass_rdt.stop();
        cond_rdt.stop();
        return 0;
//    }

    TPZLinearAnalysis an(cmesh);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    TPZSkylineStructMatrix fullstr(cmesh);
    fullstr.SetNumThreads(n_threads.get_value());
    an.SetStructuralMatrix(fullstr);
    
#ifdef PZDEBUG
    
    ofstream arg4("gmesh1.txt");
	gmesh->Print(arg4);
    
    ofstream arg3("cmesh1.txt");
	cmesh->Print(arg3);
#endif
    
    an.SetSolver(step);
    an.Run();
    total_rdt.stop();
   // an.Solution().Print("sol.txt");
    
    ResetMesh(cmesh);
    
#ifdef PZDEBUG
    ofstream arg5("cmesh_apossolve.txt");
	cmesh->Print(arg5);

    ofstream arg6("gmesh_apossolve.txt");
	gmesh->Print(arg6);
    
    ofstream file("malhageometrica.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file, true);
#endif
    
    string plotfile("Solution.vtk");
    PosProcessSol(an,plotfile);

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
		int64_t n = gmesh->NElements();
		for ( int64_t i = 0; i < n; i++ ){
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
    
    REAL diff = 1.;
	REAL conv = 0.;
	TPZVec<REAL> convdir(3,0.);
    material->SetParameters(diff, conv, convdir);
	//REAL flux = 0.;
	//material->SetInternalFlux(flux);
    
    TPZAutoPointer<TPZFunction<STATE> > force;
    force = new TPZDummyFunction<STATE>(Forcing,4);
    material->SetForcingFunction(force);
    
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
    
    TPZAutoPointer<TPZFunction<STATE> > solexata;
    solexata = new TPZDummyFunction<STATE>(SolExata,4); 
    material->SetExactSol(solexata);
    
	REAL uD=0.;
	val22(0,0)=uD;
	TPZMaterial * BCondD1 = material->CreateBC(mat1, bc2,dirichlet, val12, val22);
	cmesh->InsertMaterialObject(BCondD1);
	
	TPZMaterial * BCondD2 = material->CreateBC(mat1, bc4,dirichlet, val12, val22);
	cmesh->InsertMaterialObject(BCondD2);

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
        int64_t dim = gel->Dimension();
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
            int neq = equal.size();
            for (int eq =0; eq < neq; eq++) {
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


void Forcing(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
	double x = pt[0];
    double y = pt[1];
    disp[0]= 2.*Pi*Pi*sin(Pi*x)*sin(Pi*y);
}

void SolExata(const TPZVec<REAL> &pt, TPZVec<STATE> &disp, TPZFMatrix<STATE> &flux){
    
    disp.Resize(1, 0.);
    flux.Resize(3, 1.);
    flux(0,0)=flux(1,0)=flux(2,0)=0.;
    double x = pt[0];
    double y = pt[1];
    flux.Resize(3, 1);
    disp[0]= sin(Pi*x)*sin(Pi*y);
    flux(0,0)=-Pi*cos(Pi*x)*sin(Pi*y);
    flux(1,0)=-Pi*cos(Pi*y)*sin(Pi*x);
    flux(2,0)=2.*Pi*Pi*sin(Pi*x)*sin(Pi*y);
}

void PosProcessSol(TPZLinearAnalysis &an, std::string plotfile){
	TPZManVector<std::string,10> scalnames(2), vecnames(2);
	scalnames[0] = "Pressure";
    scalnames[1] = "ExactPressure";
    
	vecnames[0]= "Flux";
    vecnames[1]= "ExactFlux";
    
	const int dim = 2;
	int div = 0;
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an.PostProcess(div,dim);
	std::ofstream out("malha.txt");
	an.Print("nothing",out);
}



