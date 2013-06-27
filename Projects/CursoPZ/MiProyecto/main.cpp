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

int materialId = 1;
int materialBC1 = -1;

// output files  -> Because it has many energy faults
std::ofstream out("output.txt");
TPZVec<REAL> ervec(100,0.0);

void ExactSolin(const TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix<REAL> &dsol);

void BCSolin(const TPZVec<REAL> &x, TPZVec<REAL> &sol);

void Ff(const TPZVec<REAL> &x, TPZVec<REAL> &f);


void RefineGeoElements(int dim,TPZGeoMesh *gmesh,TPZManVector<REAL> &points,REAL r,REAL &distance,bool &isdefined);

TPZGeoMesh *ConstructingFicheraCorner();
TPZCompMesh *CreateMesh(TPZGeoMesh *gmesh,int dim,int hasforcingfunction);



class TCedricTest
{
public:
    
    static TPZManVector<REAL,3> fX0, fEps;
    
public:
    
    TCedricTest()
    {
        
    }

    
    void GenerateNodes(TPZGeoMesh *gmesh, int nelem);
    
    TPZGeoMesh *PyramidMesh(int nelem);
    
    TPZGeoMesh *HexaMesh(int nelem);

    void AddBoundaryElements(TPZGeoMesh *gmesh);

    TPZCompMesh *GenerateCompMesh(TPZGeoMesh *gmesh);

    static REAL fx(REAL x, REAL x0, REAL eps)
    {
        REAL result = 0.;
        REAL a = exp(-(x-x0)*(x-x0)/eps);
        REAL b = exp(-x0*x0/eps)*(1.-x);
        REAL c = exp(-(1.-x0)*(1.-x0)/eps)*x;
        result = a-b-c;
        return result;
    }

    static REAL dfx(REAL x, REAL x0, REAL eps)
    {
        REAL a = -exp(-(1-x0)*(1.-x0)/eps);
        REAL b = exp(-x0*x0/eps);
        REAL c = -2.*(x-x0)*exp(-(x-x0)*(x-x0)/eps)/eps;
        REAL result = a+b+c; 
        return result;
    }

    static REAL d2fx(REAL x, REAL x0, REAL eps)
    {
        REAL a = 2.*exp(-(x-x0)*(x-x0)/eps)/eps;
        REAL b = -4.*(x-x0)*(x-x0)*exp(-(x-x0)*(x-x0)/eps)/eps/eps;
        REAL result = a+b;
        return result;
    }
    
    static void Exact(const TPZVec<REAL> &x, TPZVec<STATE> &func, TPZFMatrix<STATE> &deriv)
    {
        REAL v[3] = {fx(x[0], fX0[0], fEps[0]),
                        fx(x[1], fX0[1], fEps[1]),
            fx(x[2], fX0[2], fEps[2])};
        func[0] = v[0]*v[1]*v[2];
        for (int i=0; i<3; i++) {
            REAL dvz = dfx(x[i], fX0[i], fEps[i]);
            deriv(i,0) = dvz*v[(i+1)%3]*v[(i+2)%3];
        }
        
    }
    
    /// verify if the faces without neighbour should be orthogonal to the main planes
    void CheckConsistency(TPZGeoMesh *mesh);
    
    void Run(int nelem)
    {
//        TPZGeoMesh *gmesh = HexaMesh(1);
        TPZGeoMesh *gmesh = PyramidMesh(nelem);
#ifdef LOG4CXX
        {
            std::stringstream sout;
            gmesh->Print(sout);
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        CheckConsistency(gmesh);
        AddBoundaryElements(gmesh);
        
        TPZCompMesh *cmesh = GenerateCompMesh(gmesh);
        
#ifdef LOG4CXX
        {
            std::stringstream sout;
            cmesh->Print(sout);
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        
        TPZAnalysis analysis(cmesh);

        analysis.SetExact(Exact);
        TPZManVector<STATE> errvec;
        analysis.PostProcessError(errvec,std::cout);

        
        TPZSkylineStructMatrix skylstr(cmesh);
        analysis.SetStructuralMatrix(skylstr);
        TPZStepSolver<STATE> step;
        step.SetDirect(ECholesky);
        analysis.SetSolver(step);
        analysis.Run();
        
        analysis.PostProcessError(errvec,std::cout);
        
        std::cout << "errvec " << errvec << std::endl;
        
    }

};

TPZManVector<REAL,3> TCedricTest::fX0(3,0.5), TCedricTest::fEps(3,0.1);

class ForceFunction : public TPZFunction<STATE>
{
    virtual void Execute(const TPZVec<STATE> &x, TPZVec<STATE> &val, TPZFMatrix<STATE> &df)
    {
        val[0] = 0.;
        for (int i=0; i<3; i++) {
            REAL vx = TCedricTest::fx(x[i], TCedricTest::fX0[i], TCedricTest::fEps[i]);
            int j = (i+1)%3;
            REAL vy = TCedricTest::fx(x[j], TCedricTest::fX0[j], TCedricTest::fEps[j]);
            int k = (j+1)%3;
            REAL vz = TCedricTest::d2fx(x[k], TCedricTest::fX0[k], TCedricTest::fEps[k]);
            val[0] -= vx*vy*vz;
        }
    }

    virtual int NFunctions()
    {
        return 1;
    }

    virtual int PolynomialOrder()
    {
        return 5;
    }

};

void formatTimeInSec(char *strtime,int timeinsec);


// MAIN FUNCTION
int main(int argc, char *argv[]) {
	
#ifdef LOG4CXX
	InitializePZLOG();
#endif
    TCedricTest cedric;
    
    cedric.Run(6);
    
    return 1;
    
    
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
	TPZMaterial *mat = new TPZMatPoisson3d(materialId,dim);
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
	TPZMaterial *bnd = mat->CreateBC(mat,materialBC1,0,val1,val2);
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

static int pyramid[2][5]=
    {
        {0,1,2,3,4},
        {4,5,6,7,2}
    };
static int tetraedra[2][4]=
{
    {1,2,5,4},
    {4,7,3,2}
};
 
void TCedricTest::GenerateNodes(TPZGeoMesh *gmesh, int nelem)
{
    gmesh->NodeVec().Resize((nelem+1)*(nelem+1)*(nelem+1));
    for (int i=0; i<=nelem; i++) {
        for (int j=0; j<=nelem; j++) {
            for (int k=0; k<=nelem; k++) {
                TPZManVector<REAL,3> x(3);
                x[0] = k*1./nelem;
                x[1] = j*1./nelem;
                x[2] = i*1./nelem;
                gmesh->NodeVec()[i*(nelem+1)*(nelem+1)+j*(nelem+1)+k].Initialize(x, *gmesh);
            }
        }
    }
}

TPZGeoMesh *TCedricTest::PyramidMesh(int nelem)
{
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    GenerateNodes(gmesh, nelem);
    
    for (int i=0; i<nelem; i++) {
        for (int j=0; j<nelem; j++) {
            for (int k=0; k<nelem; k++) {
                TPZManVector<int,8> nodes(8,0);
                nodes[0] = k*(nelem+1)*(nelem+1)+j*(nelem+1)+i;
                nodes[1] = k*(nelem+1)*(nelem+1)+j*(nelem+1)+i+1;
                nodes[2] = k*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i+1;
                nodes[3] = k*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i;
                nodes[4] = (k+1)*(nelem+1)*(nelem+1)+j*(nelem+1)+i;
                nodes[5] = (k+1)*(nelem+1)*(nelem+1)+j*(nelem+1)+i+1;
                nodes[6] = (k+1)*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i+1;
                nodes[7] = (k+1)*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i;
#ifdef LOG4CXX
                {
                    std::stringstream sout;
                    sout << "Cube nodes " << nodes;
                    LOGPZ_DEBUG(logger, sout.str())
                }
#endif
                for (int el=0; el<2; el++) 
                {
                    TPZManVector<int,5> elnodes(5);
                    for (int il=0; il<5; il++) {
                        elnodes[il] = nodes[pyramid[el][il]];
                    }
                    int index;
                    gmesh->CreateGeoElement(EPiramide, elnodes, 1, index);
                    elnodes.resize(4);
                    for (int il=0; il<4; il++) {
                        elnodes[il] = nodes[tetraedra[el][il]];
                    }
                    gmesh->CreateGeoElement(ETetraedro, elnodes, 1, index);
                }
            }
        }
    }
    gmesh->BuildConnectivity();
    return gmesh;
}

TPZGeoMesh *TCedricTest::HexaMesh(int nelem)
{
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    GenerateNodes(gmesh, nelem);
    
    for (int i=0; i<nelem; i++) {
        for (int j=0; j<nelem; j++) {
            for (int k=0; k<nelem; k++) {
                TPZManVector<int,8> nodes(8,0);
                nodes[0] = k*(nelem+1)*(nelem+1)+j*(nelem+1)+i;
                nodes[1] = k*(nelem+1)*(nelem+1)+j*(nelem+1)+i+1;
                nodes[2] = k*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i+1;
                nodes[3] = k*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i;
                nodes[4] = (k+1)*(nelem+1)*(nelem+1)+j*(nelem+1)+i;
                nodes[5] = (k+1)*(nelem+1)*(nelem+1)+j*(nelem+1)+i+1;
                nodes[6] = (k+1)*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i+1;
                nodes[7] = (k+1)*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i;
#ifdef LOG4CXX
                {
                    std::stringstream sout;
                    sout << "Cube nodes " << nodes;
                    LOGPZ_DEBUG(logger, sout.str())
                }
#endif
                int index;
                gmesh->CreateGeoElement(ECube, nodes, 1, index);
            }
        }
    }
    gmesh->BuildConnectivity();
    return gmesh;
}

/// verify if the faces without neighbour should be orthogonal to the main planes
void TCedricTest::CheckConsistency(TPZGeoMesh *mesh)
{
    int nel = mesh->NElements();
    for (int el=0; el<nel; el++) {
        TPZGeoEl *gel = mesh->ElementVec()[el];
        int nsides = gel->NSides();
        for (int is=0; is<nsides; is++) {
            TPZGeoElSide gelside(gel,is);
            if (gelside.Dimension() != 2) {
                continue;
            }
            if (gelside.Neighbour() != gelside) {
                continue;
            }
            TPZManVector<REAL,2> xi(2,0.);
            gelside.CenterPoint(xi);
            TPZFNMatrix<6,REAL> axes(2,3);
            TPZFNMatrix<4,REAL> jac(2,2),jacinv(2,2);
            REAL detjac;
            gelside.Jacobian(xi, jac, axes, detjac, jacinv);
            TPZManVector<REAL,3> x(3,0.);
            gelside.X(xi, x);
            TPZManVector<REAL,3> normal(3);
            normal[0] = fabs(axes(0,1)*axes(1,2)-axes(0,2)*axes(1,1));
            normal[1] = fabs(-axes(0,0)*axes(1,2)+axes(0,2)*axes(1,0));
            normal[2] = fabs(axes(0,0)*axes(1,1)-axes(0,1)*axes(1,0));
            REAL tol = 1.e-6;
            REAL xmin = 1., xmax = 0.;
            int numtol = 0;
            for (int i=0; i<3; i++) {
                if(xmin > x[i]) xmin = x[i];
                if (xmax < x[i]) {
                    xmax = x[i];
                }
                if (normal[i] > tol) {
                    numtol++;
                }
            }
            if (numtol != 1) {
                DebugStop();
            }
            if (xmin > tol && xmax < 1.-tol) {
                DebugStop();
            }
        }
    }
}

TPZCompMesh *TCedricTest::GenerateCompMesh(TPZGeoMesh *gmesh)
{
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    int dim=3;
    
    TPZMatPoisson3d *poiss = new TPZMatPoisson3d(1,dim);
    
    TPZAutoPointer<TPZFunction<STATE> > force = new ForceFunction;
    
    poiss->SetForcingFunction(force );
    
    TPZFMatrix<STATE> val1(1,1,0.), val2(1,1,0.);
    TPZBndCond *bc = new TPZBndCond(poiss, -1, 0, val1, val2);
    
    cmesh->InsertMaterialObject(poiss);
    cmesh->InsertMaterialObject(bc);
    
    cmesh->AutoBuild();
    
    return cmesh;
}

void TCedricTest::AddBoundaryElements(TPZGeoMesh *gmesh)
{
    int nelem = gmesh->NElements();
    for (int el = 0; el<nelem; el++) {
        TPZGeoEl *gel = gmesh->ElementVec()[el];
        int ns = gel->NSides();
        for (int is=0; is<ns; is++) {
            TPZGeoElSide gelside(gel,is);
            if (gelside.Dimension() != 2) {
                continue;
            }
            if (gelside.Neighbour() == gelside) {
                TPZGeoElBC(gelside, -1);
            }
        }
    }
}

