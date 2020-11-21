
#include "TPZGenGrid2D.h"
#include "TPZGenGrid3D.h"
#include "pzcheckgeom.h"
#include "pzgmesh.h"

#include "TPZPrimalPoisson.h"
#include "TPZMatLaplacianHybrid.h"
#include "mixedpoisson.h"

#include "pzbndcond.h"
#include "TPZAnalyticSolution.h"
#include "pzanalysis.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"

#include "TPZMHMeshControl.h"
#include "TPZMHMixedMeshControl.h"
#include "TPZMHMixedHybridMeshControl.h"
#include "pzlog.h"


#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.testmhm"));
#endif

#ifdef USING_BOOST

#ifndef WIN32
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MAIN pz meshHDiv tests

#include <boost/test/unit_test.hpp>

// generate the geometric mesh
TPZGeoMesh *GenerateGeoMesh(int dim, int numel, int ndiv);

// generate the MHM H1 mesh
TPZMHMeshControl *GenerateMHMH1(int dim, int numel, int ndiv);

// generate the MHM HDiv mesh
TPZMHMixedMeshControl *GenerateMHMHDiv(int dim, int numel, int ndiv);

// generate the MHM Hybrid HDiv mesh
TPZMHMixedHybridMeshControl *GenerateMHMHDivHybrid(int dim, int numel, int ndiv);

// solve the problem and verify the errors
bool SolveProblem(TPZCompMesh *cmesh);

// insert H1 material in the multiphysics mesh
void InsertH1Materials(TPZAutoPointer<TPZMultiphysicsCompMesh> cmesh);

// insert HDiv material in the multiphysics mesh
void InsertHDivMaterials(TPZAutoPointer<TPZMultiphysicsCompMesh> cmesh);

#ifdef _AUTODIFF
TLaplaceExample1 LaplaceExact;
#endif

// Tests for the 'voidflux' class.
BOOST_AUTO_TEST_SUITE(mmh_tests)

BOOST_AUTO_TEST_CASE(MHMHDivHybrid)
{
    #ifdef LOG4CXX
        InitializePZLOG();
    #endif

    std::cout << "Running MHM Hybrid HDiv test\n";
    int dim = 2;
    int numelcoarse = 3;
    int ndiv = 1;
    TPZMHMixedHybridMeshControl *mhm = GenerateMHMHDivHybrid(dim, numelcoarse, ndiv);
    TPZCompMesh *cmesh = mhm->CMesh().operator->();
    bool correct = SolveProblem(cmesh);
    delete mhm;
    
    BOOST_CHECK(correct);
}



BOOST_AUTO_TEST_CASE(MHMH1)
{
    std::cout << "Running MHM H1 test\n";
    TPZMHMeshControl *mhm = GenerateMHMH1(2, 3, 1);
    TPZCompMesh *cmesh = mhm->CMesh().operator->();
    bool correct = SolveProblem(cmesh);
    delete mhm;
    
    BOOST_CHECK(correct);
}


BOOST_AUTO_TEST_CASE(MHMHDiv)
{

    std::cout << "Running MHM HDiv test\n";
    TPZMHMixedMeshControl *mhm = GenerateMHMHDiv(2, 3, 1);
    TPZCompMesh *cmesh = mhm->CMesh().operator->();
    bool correct = SolveProblem(cmesh);
    delete mhm;
    
    BOOST_CHECK(correct);
}


BOOST_AUTO_TEST_SUITE_END()

#endif

// solve the problem and verify the errors
bool SolveProblem(TPZCompMesh *cmesh)
{
    bool renumber = true;
    TPZAnalysis an(cmesh,renumber);
    TPZSkylineStructMatrix str(cmesh);
    an.SetStructuralMatrix(str);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);
    LaplaceExact.fExact = TLaplaceExample1::EX;
    an.SetExact(LaplaceExact.ExactSolution());
    cmesh->MaterialVec()[1]->SetForcingFunction(LaplaceExact.ForcingFunction());
    cmesh->MaterialVec()[-1]->SetForcingFunction(LaplaceExact.Exact());
//    an.Solution().Print("Sol =",std::cout, EMathematicaInput);
    an.Run();
    
#ifdef PZDEBUG
    TPZStack<std::string> scalnames,vecnames;
    scalnames.Push("Pressure");
    std::string outfile("MHMPlot.vtk");
    an.DefineGraphMesh(cmesh->Dimension(), scalnames, vecnames, outfile);
    an.PostProcess(1);
//    an.Solution().Print("Solution");
#endif
    
    TPZVec<REAL> errors(3,0.);
    int64_t nel = cmesh->NElements();
    cmesh->ElementSolution().Resize(nel, 3);
    an.PostProcessError(errors,false);
    std::cout << "Errors " << errors << std::endl;
    bool correct = (errors[2] < 1.e-6);
    return correct;
}



// generate the gometric mesh
TPZGeoMesh *GenerateGeoMesh(int dim, int numel, int ndiv)
{
//    TPZGenGrid2D(const TPZVec<int> &nx, const TPZVec<REAL> &x0, const TPZVec<REAL> &x1, int numl = 1, REAL rot = 0.5);
    TPZGeoMesh *gmesh = 0;
    if(dim == 2)
    {
        gmesh = new TPZGeoMesh();
        gmesh->SetDimension(2);
        TPZVec<int> nx = {numel,numel};
        TPZVec<REAL> x0 = {0.,0.,0.}, x1 = {1.,1.,0.};
        TPZGenGrid2D gengrid(nx,x0,x1);
        gengrid.Read(gmesh);
        gengrid.SetBC(gmesh, 4, -1);
        gengrid.SetBC(gmesh, 5, -1);
        gengrid.SetBC(gmesh, 6, -1);
        gengrid.SetBC(gmesh, 7, -1);
    }
    else if(dim == 3)
    {
        TPZVec<REAL> minX = {0.,0.,0.}, maxX = {1.,1.,1.};
        TPZVec<int> nelDiv = {numel,numel,numel};
        MMeshType eltype = MMeshType::EHexahedral;
        TPZGenGrid3D gengrid(minX, maxX, nelDiv, eltype);
        gengrid.BuildVolumetricElements(1);
        gmesh = gengrid.BuildBoundaryElements(-1, -1, -1, -1, -1, -1);
    }
    
    TPZCheckGeom check(gmesh);
    check.UniformRefine(ndiv);
    return gmesh;
}

void PropagatePartition(TPZGeoEl *gel, TPZVec<int64_t> &partitionindex)
{
    int64_t index = gel->Index();
    int64_t partition = partitionindex[index];
    if(partition == -1) DebugStop();
    TPZManVector<TPZGeoEl *> subels;
    int HasSubEl = gel->HasSubElement();
    if(HasSubEl)
    {
        int nsubel = gel->NSubElements();
        for (int sub=0; sub<nsubel; sub++) {
            TPZGeoEl *subel =gel->SubElement(sub);
            int64_t subelindex = subel->Index();
            partitionindex[subelindex] = partition;
            PropagatePartition(subel, partitionindex);
        }
    }
}

// generate the MHM H1 mesh
TPZMHMeshControl *GenerateMHMH1(int dim, int numel, int ndiv)
{
    TPZGeoMesh *gmesh = GenerateGeoMesh(dim, numel, ndiv);
    TPZMHMeshControl *mhm = new TPZMHMeshControl(gmesh);
    int64_t nel = gmesh->NElements();
    int64_t nelcoarse = numel*numel;
    if(dim == 3) nelcoarse *=numel;
    TPZVec<int64_t> partition(nelcoarse,-1);
    for(int el=0; el<nelcoarse; el++)
    {
        partition[el] = el;
    }
    std::cout << "Partition " << partition << std::endl;
    mhm->DefinePartitionbyCoarseIndices(partition);
    InsertH1Materials(mhm->CMesh());
    mhm->SetInternalPOrder(2);
    mhm->SetSkeletonPOrder(1);
    mhm->fMaterialIds.insert(1);
    mhm->fMaterialBCIds.insert(-1);
    mhm->SetLagrangeAveragePressure(true);
    mhm->BuildComputationalMesh(true);
    return mhm;
}

// generate the MHM HDiv mesh
TPZMHMixedMeshControl *GenerateMHMHDiv(int dim, int numel, int ndiv)
{
    TPZGeoMesh *gmesh = GenerateGeoMesh(dim, numel, ndiv);
    TPZMHMixedMeshControl *mhm = new TPZMHMixedMeshControl(gmesh);
    int64_t nel = gmesh->NElements();
    int64_t nelcoarse = numel*numel;
    if(dim == 3) nelcoarse *=numel;
    TPZVec<int64_t> partition(nelcoarse,-1);
    for(int el=0; el<nelcoarse; el++)
    {
        partition[el] = el;
    }
    std::cout << "Partition " << partition << std::endl;
    mhm->DefinePartitionbyCoarseIndices(partition);
    InsertHDivMaterials(mhm->CMesh());
    mhm->SetInternalPOrder(2);
    mhm->SetSkeletonPOrder(1);
    mhm->fMaterialIds.insert(1);
    mhm->fMaterialBCIds.insert(-1);
    mhm->SetLagrangeAveragePressure(true);
    mhm->BuildComputationalMesh(true);
    return mhm;

}

// generate the MHM Hybrid HDiv mesh
TPZMHMixedHybridMeshControl *GenerateMHMHDivHybrid(int dim, int numel, int ndiv)
{
    TPZGeoMesh *gmesh = GenerateGeoMesh(dim, numel, ndiv);
    TPZMHMixedHybridMeshControl *mhm = new TPZMHMixedHybridMeshControl(gmesh);
    int64_t nel = gmesh->NElements();
    int64_t nelcoarse = numel*numel;
    if(dim == 3) nelcoarse *=numel;
    TPZVec<int64_t> partition(nelcoarse,-1);
    for(int el=0; el<nelcoarse; el++)
    {
        partition[el] = el;
    }
    std::cout << "Partition " << partition << std::endl;
    mhm->DefinePartitionbyCoarseIndices(partition);
    InsertHDivMaterials(mhm->CMesh());
    mhm->SetInternalPOrder(1);
    mhm->SetSkeletonPOrder(1);
    mhm->fMaterialIds.insert(1);
    mhm->fMaterialBCIds.insert(-1);
    mhm->SetLagrangeAveragePressure(true);
    bool substructure = true;
    mhm->BuildComputationalMesh(substructure);
    return mhm;

}




// insert H1 material in the multiphysics mesh
void InsertH1Materials(TPZAutoPointer<TPZMultiphysicsCompMesh> cmesh)
{
    TPZMatLaplacianHybrid *mat = new TPZMatLaplacianHybrid(1,cmesh->Dimension());
//    TPZBndCond(TPZMaterial * material,int id,int type,const TPZFMatrix<STATE> &val1,const TPZFMatrix<STATE> &val2) :
    TPZFNMatrix<2> val1(1,1,0.),val2(1,1,1.);
    TPZBndCond *bc = new TPZBndCond(mat,-1,0,val1,val2);
    cmesh->InsertMaterialObject(mat);
    cmesh->InsertMaterialObject(bc);
}

// insert H1 material in the multiphysics mesh
void InsertHDivMaterials(TPZAutoPointer<TPZMultiphysicsCompMesh> cmesh)
{
    TPZMixedPoisson *mat = new TPZMixedPoisson(1, cmesh->Dimension()); //Using standard
//    TPZBndCond(TPZMaterial * material,int id,int type,const TPZFMatrix<STATE> &val1,const TPZFMatrix<STATE> &val2) :
    TPZFNMatrix<2> val1(1,1,0.),val2(1,1,1.);
    TPZBndCond *bc = new TPZBndCond(mat,-1,0,val1,val2);
    cmesh->InsertMaterialObject(mat);
    cmesh->InsertMaterialObject(bc);
}

