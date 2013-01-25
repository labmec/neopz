/**
 * @file
 * @brief Projeto elaborado para encontrar os errores seguintes:
 * - Cuando refina elementos tridimensionais existem fDepend dos connects que não foram previamente deletados
 * - Para uma ordem alta em 3D, no TPZAnalysis, posprocessing pega dimensão maior do Block Information
 * - No pzintel.h acontece o DebugStop no check que o Philippe introduziu para pegar problemas com a ordem de interpolação
 */

#include "pzshapelinear.h"

#include "pzgengrid.h"
#include "TPZExtendGridDimension.h"
#include "TPZVTKGeoMesh.h"

#include "pzbstrmatrix.h"

#include "pzintel.h"
#include "pzcompel.h"
#include "pzcheckmesh.h"

#include "pzmaterial.h"
#include "pzbndcond.h"
#include "pzelasmat.h"
#include "pzpoisson3d.h"

#include "pzanalysis.h"
#include "pzstepsolver.h"

#include "TPZRefPatternTools.h"

#include "TPZParSkylineStructMatrix.h"

#include <stdio.h>

#include "pzlog.h"

#include "pzgeoelbc.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.fischera"));
#endif

using namespace std;
using namespace pzshape;


// output files  -> Because it has many energy faults
std::ofstream out("output.txt");
TPZVec<REAL> ervec(100,0.0);

void ExactSolin(const TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix<REAL> &dsol);

void BCSolin(const TPZVec<REAL> &x, TPZVec<REAL> &sol);

void Ff(const TPZVec<REAL> &x, TPZVec<REAL> &f);


void RefineGeoElements(int dim,TPZGeoMesh *gmesh,TPZManVector<REAL> &points,REAL r,REAL &distance,bool &isdefined);

TPZGeoMesh *ConstructingFicheraCorner();
TPZCompMesh *CreateMesh(TPZGeoMesh *gmesh,int dim,int hasforcingfunction);

void formatTimeInSec(char *strtime,int timeinsec);


// MAIN FUNCTION
int main(int argc, char *argv[]) {
	
#ifdef LOG4CXX
	InitializePZLOG();
#endif
	
	// Initializing a ref patterns
//	gRefDBase.InitializeAllUniformRefPatterns();
	gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
	gRefDBase.InitializeUniformRefPattern(ECube);
    

	//-----------  INITIALIZING CONSTRUCTION OF THE MESHES
	int i, nref, NRefs = 5;
	int dim = 3;
	int nelem = 0;
    
    for(nref=2;nref<NRefs;nref++) {
		
        // Constructing geometric mesh as Fichera corner using hexahedra
        TPZGeoMesh *gmesh3D = ConstructingFicheraCorner();
// h_refinement
        TPZManVector<REAL> point(3,0.);
        REAL r = 0.0, radius = 0.9;
        bool isdefined = false;
 
		// The initial mesh is a cube then we are going to refine one level (8 cubes will be created)
		RefineGeoElements(3,gmesh3D,point,r,radius,isdefined);
		TPZVec<TPZGeoEl*> sub;
		// The fifth element will be divided
		gmesh3D->ElementVec()[4]->Divide(sub);
        
        // Creating computational mesh
        /** Set polynomial order */
        int p = 2;
        TPZCompEl::SetgOrder(p);
        TPZCompMesh *cmesh = CreateMesh(gmesh3D,dim,1);
		cmesh->SetName("Computational mesh for Fichera problem");
		dim = cmesh->Dimension();
		cmesh->Print();

// p-refinement
        TPZInterpolatedElement *intel;
		nelem = cmesh->NElements();
		// Searching a first cube not null and the order will be decremented
        for(i=0;i<nelem;i++) {
            intel = (TPZInterpolatedElement*)(cmesh->ElementVec()[i]);
            if(!intel) continue;
			intel->PRefine(1);
			break;
        }
        cmesh->ExpandSolution();
        cmesh->CleanUpUnconnectedNodes();

		//--- END construction of the meshes
        
		/** Variable names for post processing */
        TPZStack<std::string> scalarnames, vecnames;
		scalarnames.Push("POrder");
		scalarnames.Push("Solution");
        
        // Introduzing exact solution
        TPZAnalysis an (cmesh);
        TPZSkylineStructMatrix strskyl(cmesh);
        an.SetStructuralMatrix(strskyl);
        
        TPZStepSolver<REAL> *direct = new TPZStepSolver<REAL>;
        direct->SetDirect(ECholesky);
        an.SetSolver(*direct);
        delete direct;
        direct = 0;
        
        // Solving
        an.Run();
               
        // Post processing
        std::string filename = "Poisson3DSol.vtk";
        an.DefineGraphMesh(dim,scalarnames,vecnames,filename);
        
        an.PostProcess(0,dim);
        
        delete cmesh;
        delete gmesh3D;
	}
    out.close();
	return 0;
}

////////////////////////////////////////////////////////////////////////////////////////
//////////   FICHERA CORNER - Problem as Anders Solin Presentation   ///////////////////
////////////////////////////////////////////////////////////////////////////////////////

TPZGeoMesh *ConstructingFicheraCorner() {
    REAL co[8][3] = {
        {0.,0.,0.},
        {1.,0.,0.},
        {1.,1.,0.},
        {0.,1.,0.},
        {0.,0.,1.},
        {1.,0.,1.},
        {1.,1.,1.},
        {0.,1.,1.}
    };
    int indices[1][8] = {{0,1,2,3,4,5,6,7}};
    
    const int nelem = 1;
    int nnode = 8;
    
    TPZGeoEl *elvec[nelem];
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    
    int nod;
    for(nod=0; nod<nnode; nod++) {
        int nodind = gmesh->NodeVec().AllocateNewElement();
        TPZVec<REAL> coord(3);
        coord[0] = co[nod][0];
        coord[1] = co[nod][1];
        coord[2] = co[nod][2];
        gmesh->NodeVec()[nodind] = TPZGeoNode(nod,coord,*gmesh);
    }
    
    int el;
    for(el=0; el<nelem; el++) {
        TPZManVector<int> nodind(8);
        for(nod=0; nod<8; nod++) nodind[nod]=indices[el][nod];
        int index;
        elvec[el] = gmesh->CreateGeoElement(ECube,nodind,1,index);
    }
    
    gmesh->BuildConnectivity();
    
 //   TPZGeoElBC gbc1(elvec[0],20,-1);
 //   TPZGeoElBC gbc2(elvec[0],25,-2);
    return gmesh;
}

void ExactSolin(const TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix<REAL> &dsol) {
	REAL quad_r = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
	sol[0] = sqrt( sqrt (quad_r) );
	if(!IsZero(sol[0])) {
		REAL den = sol[0]*sol[0]*sol[0];
		dsol(0,0) = .5*x[0]/den;
		dsol(1,0) = .5*x[1]/den;
		dsol(1,0) = .5*x[2]/den;
	}
	else {
		dsol(0,0) = dsol(1,0) = dsol(2,0) = 0.;
	}
}

void BCSolin(const TPZVec<REAL> &x, TPZVec<REAL> &bcsol) {
	REAL quad_r = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
	bcsol[0] = sqrt( sqrt (quad_r) );	
}

void Ff(const TPZVec<REAL> &x, TPZVec<REAL> &f) {
	REAL quad_r = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
	REAL raiz = sqrt( sqrt(quad_r));
	f[0] = -3./(4.*(raiz*raiz*raiz));
}

TPZCompMesh *CreateMesh(TPZGeoMesh *gmesh,int dim,int hasforcingfunction) {
	
	TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
	cmesh->SetDefaultOrder(TPZCompEl::GetgOrder());
	cmesh->SetAllCreateFunctionsContinuous();
	
	// Creating Poisson material
	TPZMaterial *mat = new TPZMatPoisson3d(1,dim);
	TPZVec<REAL> convd(3,0.);
	((TPZMatPoisson3d *)mat)->SetParameters(1.,0.,convd);
	if(hasforcingfunction) {
		mat->SetForcingFunction(new TPZDummyFunction<STATE>(Ff));
	}
	cmesh->InsertMaterialObject(mat);
	// Make compatible dimension of the model and the computational mesh
	cmesh->SetDimModel(mat->Dimension());
	cmesh->SetAllCreateFunctionsContinuous();
    
	// Boundary conditions
	// Dirichlet
	TPZAutoPointer<TPZFunction<STATE> > FunctionBC = new TPZDummyFunction<STATE>(BCSolin);
	TPZFMatrix<REAL> val1(dim,dim,0.),val2(dim,1,0.);
	TPZMaterial *bnd = mat->CreateBC(mat,-1,0,val1,val2);
	bnd->SetForcingFunction(FunctionBC);
	cmesh->InsertMaterialObject(bnd);
	
	cmesh->AutoBuild();
    
#ifdef LOG4CXX
    {
        std::stringstream sout;
        TPZCheckMesh tst(cmesh,&sout);
        tst.VerifyAllConnects();
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif    
    
    cmesh->AdjustBoundaryElements();
    cmesh->ExpandSolution();
	cmesh->CleanUpUnconnectedNodes();
	return cmesh;
}

void RefineGeoElements(int dim,TPZGeoMesh *gmesh,TPZManVector<REAL> &point,REAL r,REAL &distance,bool &isdefined) {
	TPZManVector<REAL> centerpsi(3), center(3);
	// Refinamento de elementos selecionados
	TPZGeoEl *gel;
	TPZVec<TPZGeoEl *> sub;
	
	int nelem = 0;
	int ngelem=gmesh->NElements();
	// na esquina inferior esquerda Nó = (0,-1,0)
	while(nelem<ngelem) {
		gel = gmesh->ElementVec()[nelem++];
		if(gel->Dimension()!=dim || gel->HasSubElement()) continue;
		gel->CenterPoint(gel->NSides()-1,centerpsi);
		gel->X(centerpsi,center);
		if(!isdefined) {
			TPZVec<REAL> FirstNode(3,0.);
			gel->CenterPoint(0,centerpsi);
			gel->X(centerpsi,FirstNode);
			REAL distancia = TPZGeoEl::Distance(center,FirstNode);
			if(distancia > distance) distance = distancia;
			isdefined = true;
		}
		REAL centerdist = TPZGeoEl::Distance(center,point);
		if(fabs(r-centerdist) < distance) {
			gel->Divide(sub);
		}
	}
}

