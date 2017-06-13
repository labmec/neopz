
#include "pzlog.h"
#include "pzgmesh.h"
#include "pzmanvector.h"
#include "TPZMHMeshControl.h"
#include "TPZVTKGeoMesh.h"

#include "pzelasmat.h"
#include "TPZElasticity2DHybrid.h"
#include "pzmat1dlin.h"
#include "pzbndcond.h"
#include "pzanalysis.h"

#include "meshgen.h"

/// Insert material objects for the MHM Mesh solution
void InsertMaterialObjects(TPZMHMeshControl &control);

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mainskeleton"));
#endif


const int matInterno = 1;
const int matCoarse = 2;
const int skeleton = 4;
const int secondskeleton = 3;
const int matInterface = 5;
const int matpressure = 6;

const int dirichlet = 0;
const int neumann = 1;
const int mixed = 2;

int const bc1=-1;
int const bc2=-2;
int const bc3=-3;
int const bc4=-4;
int const bc5=-5;

struct TRunConfig
{
    int nelxref;
    int nelyref;
    int numHDivisions;
    int pOrder;
};



int main(int argc, char *argv[])
{
    TExceptionManager except;

    /// computation type :
    // (0) - compute reference mesh
    // (1) - compute MHM H1 mesh and compute MHM(div) mesh
    int ComputationType = 1;
    /// numhdiv - number of h-refinements
    int NumHDivision = 1;
    /// PolynomialOrder - p-order
    int PolynomialOrder = 2;
    
    TRunConfig Configuration;
    
    TPZManVector<REAL,3> divsigma(2,0.), x(2,0.5);
    TElasticityExample1::DivSigma(x, divsigma);
    std::cout << "x = " << x << " divsigma = " << divsigma << std::endl;
    
    Configuration.pOrder = PolynomialOrder;
    Configuration.numHDivisions = NumHDivision;

    if (argc == 4)
    {
        ComputationType = atoi(argv[1]);
        NumHDivision = atoi(argv[2]);
        PolynomialOrder = atoi(argv[3]);
    }
    HDivPiola = 1;
#ifdef LOG4CXX
    InitializePZLOG();
#endif

    // verifying differences between the MHM-original and MHM with mixed approximations
    int nelx = 1;
    int nely = 1;
    {
        std::ofstream out("DiffResults.nb",std::ios::app);
        out << "(* Running quadrilateral elastic mesh with numsubdomains " << nelx << ", " << nely << " *)\n";
    }
    TPZGeoMesh *gmesh = 0;
    TPZVec<long> coarseindices;
   
    TPZManVector<REAL,3> x0(3,0.),x1(3,1.);
    x1[2] = 0.;
    int ndiv = NumHDivision;
    gmesh = MalhaGeomFredQuadrada(nelx, nely, x0, x1, coarseindices, ndiv);

    TPZAutoPointer<TPZGeoMesh> gmeshauto(gmesh);
    TPZAutoPointer<TPZMHMeshControl> MHM;
    
    
    {
        TPZAutoPointer<TPZGeoMesh> gmeshauto = new TPZGeoMesh(*gmesh);
        TPZMHMeshControl *mhm = new TPZMHMeshControl(gmeshauto,coarseindices);
        MHM = mhm;
        TPZMHMeshControl &meshcontrol = *mhm;
        
        meshcontrol.SetLagrangeAveragePressure(false);
        
        InsertMaterialObjects(meshcontrol);
        
        int porder = PolynomialOrder;
        // to avoid singular internal matrices
        if (NumHDivision == 0 && porder == 1) {
            porder++;
        }
        meshcontrol.SetInternalPOrder(porder);
        meshcontrol.SetSkeletonPOrder(2);
        
        meshcontrol.CreateSkeletonElements(skeleton);
        
        meshcontrol.DivideSkeletonElements(0);
        //        meshcontrol.Hybridize(secondskeleton, matpressure);
        
        bool substructure = false;
        meshcontrol.BuildComputationalMesh(substructure);
#ifdef PZDEBUG
        if(1)
        {
            std::ofstream file("GMeshControl.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(meshcontrol.GMesh().operator->(), file);
        }
#endif
#ifdef PZDEBUG
        if(1)
        {
            std::ofstream out("MHMMeshControl.txt");
            meshcontrol.Print(out);
        }
#endif
        
        std::cout << "MHM Computational meshes created\n";
#ifdef PZDEBUG
        if(1)
        {
            std::ofstream gfile("geometry.txt");
            gmesh->Print(gfile);
            
            std::ofstream out_mhm("MHM_h1.txt");
            meshcontrol.CMesh()->Print(out_mhm);
            
        }
#endif
        std::cout << "Number of equations MHM equals " << MHM->CMesh()->NEquations() << std::endl;
        
    }
    std::string configuration;
    
    {
        std::stringstream sout;
        sout << "H" << NumHDivision << "-P" << PolynomialOrder;
        configuration = sout.str();
    }

    // compute the MHM solution
    SolveProblem(MHM->CMesh(), MHM->GetMeshes(), "MHMElast", configuration);


    TPZAnalysis locanalysis(MHM->CMesh(),false);
    locanalysis.SetExact(TElasticityExample1::GradU);
    TPZVec<STATE> errors(3,0.);
    locanalysis.PostProcessError(errors);
    std::cout << "Errors computed " << errors << std::endl;
    return 0;
}

void InsertMaterialObjects(TPZMHMeshControl &control)
{
    TPZCompMesh &cmesh = control.CMesh();
    /// criar materiais
//    int dim = cmesh.Dimension();
    control.SetProblemType(TPZMHMeshControl::EElasticity2D);
    STATE Young = 1000., nu = 0.3, fx = 0., fy = 0.;
    TPZElasticity2DHybrid *material1 = new TPZElasticity2DHybrid(matInterno,Young,nu,fx,fy);
    material1->SetPlaneStrain();
    
    material1->SetForcingFunction(TElasticityExample1::ForcingFunction());
    TPZMaterial * mat1(material1);
    
    TPZMat1dLin *materialCoarse = new TPZMat1dLin(matCoarse);
    TPZFNMatrix<1,STATE> xk(2,2,0.),xb(2,2,0.),xc(2,2,0.),xf(2,1,0.);
    materialCoarse->SetMaterial(xk, xc, xb, xf);
    
    cmesh.InsertMaterialObject(materialCoarse);
    materialCoarse = new TPZMat1dLin(skeleton);
    materialCoarse->SetMaterial(xk, xc, xb, xf);
    cmesh.InsertMaterialObject(materialCoarse);
    materialCoarse = new TPZMat1dLin(secondskeleton);
    materialCoarse->SetMaterial(xk, xc, xb, xf);
    cmesh.InsertMaterialObject(materialCoarse);
    materialCoarse = new TPZMat1dLin(matpressure);
    materialCoarse->SetMaterial(xk, xc, xb, xf);
    cmesh.InsertMaterialObject(materialCoarse);
    
    
    
    //    REAL diff = -1.;
    //	REAL conv = 0.;
    //	TPZVec<REAL> convdir(3,0.);
    //	REAL flux = 8.;
    //
    //	material1->SetParameters(diff, conv, convdir);
    //	material1->SetInternalFlux( flux);
    
    cmesh.InsertMaterialObject(mat1);
    
    
    control.SwitchLagrangeMultiplierSign(true);
    ///Inserir condicao de contorno
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
    //BC -1
    val1(0,0) = 0.;
    val2.Zero();
    val1(0,0) = 0;
    val1(1,1) = 0;
    TPZMaterial * BCondD1 = material1->CreateBC(mat1, bc1,dirichlet, val1, val2);
    BCondD1->SetForcingFunction(TElasticityExample1::DirichletFunction());
    cmesh.InsertMaterialObject(BCondD1);
    
    //BC -2
    val1.Zero();
    val2(0,0) = 10.;
    TPZMaterial * BCondD2 = material1->CreateBC(mat1, bc2,dirichlet, val1, val2);
    BCondD2->SetForcingFunction(TElasticityExample1::DirichletFunction());
    cmesh.InsertMaterialObject(BCondD2);
    
    //BC -3
    val1.Zero();
    val2.Zero();
    TPZMaterial * BCondD3 = material1->CreateBC(mat1, bc3,dirichlet, val1, val2);
    BCondD3->SetForcingFunction(TElasticityExample1::DirichletFunction());
    cmesh.InsertMaterialObject(BCondD3);
    
    //BC -4
    val1(0,0) = 0;
    val1(1,1) = 1.e9;
    val2(0,0) = -1.;
    TPZMaterial * BCondD4 = material1->CreateBC(mat1, bc4,dirichlet, val1, val2);
    BCondD4->SetForcingFunction(TElasticityExample1::DirichletFunction());
    cmesh.InsertMaterialObject(BCondD4);
    
    //BC -5: dirichlet nulo
    TPZMaterial * BCondD5 = material1->CreateBC(mat1, bc5,dirichlet, val1, val2);
    cmesh.InsertMaterialObject(BCondD5);
}
