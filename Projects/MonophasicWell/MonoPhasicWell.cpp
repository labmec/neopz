
#include "tpzautopointer.h"
#include "pzlog.h"
#include <time.h>

#include "tpzgeoelrefpattern.h"
#include "TPZRefPatternTools.h"
#include "TPZGeoLinear.h"
#include "pzgeopoint.h"

#include "tpzhierarquicalgrid.h"
#include "TPZVTKGeoMesh.h"
#include "TPZMonoPhaseWell.h"
#include "pzfstrmatrix.h"

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzbndcond.h"

#include "TPZVTKGeoMesh.h"
#include "TPZReadGIDGrid.h"

#include "pznonlinanalysis.h"
#include "pzstepsolver.h"
#include "TPZSkylineNSymStructMatrix.h"



#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.WellFlow"));
#endif

void ParametricfunctionS(const TPZVec<REAL> &par, TPZVec<REAL> &X);
void Ffunction(const TPZVec<REAL> &pt, TPZVec<STATE> &ff);

TPZGeoMesh * WellMesh(REAL s, REAL ds, int nelements);
void PrintGeoMesh(TPZGeoMesh * gmesh);

TPZCompMesh * CmeshFlux(int qorder, TPZGeoMesh * gmesh);
TPZCompMesh * CmeshPressure(int porder, TPZGeoMesh * gmesh);
TPZCompMesh * CmeshMixed(TPZGeoMesh * gmesh);
TPZCompMesh * CmeshMixedInitial(TPZGeoMesh * gmesh);
void PostProcesswithVTK(TPZAnalysis *an);
void PrintLS(TPZAnalysis *an);
void AssembleLastState(TPZAnalysis *an, TPZCompMesh * cmesh);
void AssembleNextState(TPZAnalysis *an, TPZCompMesh * cmesh);
void NewtonIterations(TPZAnalysis *an, TPZManVector<TPZCompMesh *> meshvector, TPZCompMesh * cmesh);
void TimeForward(TPZAnalysis *an, TPZManVector<TPZCompMesh *> meshvector, TPZCompMesh * cmesh);

void InitialProblem(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh * > meshvec);
void TransientProblem(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh * > meshvec);

int main()
{
    std::string dirname = PZSOURCEDIR;
    gRefDBase.InitializeUniformRefPattern(EOned);
    
#ifdef PZDEBUG
#ifdef LOG4CXX
    
    std::string FileName = dirname;
    FileName = dirname + "/Projects/MonophasicWell/";
    FileName += "WellFlowLog.cfg";
    InitializePZLOG(FileName);
    
#endif
#endif
    
    // Geometry of well
    REAL s= 0.0;
    REAL ds = 10.0;
    int nelements = 20;
    TPZGeoMesh * gmesh = WellMesh(s, ds, nelements);
    PrintGeoMesh(gmesh);
    
    // Computational Mesh
    int q = 2;
    int p = 1;
    
    TPZManVector<TPZCompMesh *>meshvec(2);
    meshvec[0] = CmeshFlux(q,gmesh);
    meshvec[1] = CmeshPressure(p,gmesh);
    
    InitialProblem(gmesh, meshvec);
    
    TransientProblem(gmesh, meshvec);
    
    std::cout << " Process complete normally." << std::endl;
    return 0;
}

void InitialProblem(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh * > meshvec)
{
    TPZCompMesh * cmeshwell = CmeshMixedInitial(gmesh);
    
    int wellmatid = 1;
    REAL dt = 1.0e14;
    TPZMaterial * Matwell = cmeshwell->FindMaterial(wellmatid);
    TPZMonoPhaseWell *Well = dynamic_cast<TPZMonoPhaseWell *>(Matwell);
    Well->Setdt(dt);
    
    // Transferindo para a multifisica
    TPZBuildMultiphysicsMesh::AddElements(meshvec, cmeshwell);
    TPZBuildMultiphysicsMesh::AddConnects(meshvec, cmeshwell);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, cmeshwell);
    
#ifdef PZDEBUG
    std::ofstream dumpfile("ComputationaMeshMultiphysicInitial.txt");
    cmeshwell->Print(dumpfile);
#endif
    
    
    // Setting up the analysis configuration
    int numofThreads = 4;
    TPZAnalysis * an = new TPZAnalysis(cmeshwell,true);
    TPZSkylineNSymStructMatrix skylnsym(cmeshwell);
    TPZStepSolver<STATE> step;
    skylnsym.SetNumThreads(numofThreads);
    step.SetDirect(ELU);
    an->SetStructuralMatrix(skylnsym);
    an->SetSolver(step);
    
    int64_t numofequ = an->Solution().Rows();
    std::cout << " Number of DOF " << numofequ << std::endl;
    
    // Initialize solution with hydrostatic gradient
    std::cout << "Time Value (days):  " << 0.0 << std::endl;
    std::cout<<  "Time step: " << 0 << std::endl;
    
    const clock_t tinia = clock();
    NewtonIterations(an,meshvec,cmeshwell);
    const clock_t tenda = clock();
    const REAL timea = REAL(REAL(tenda - tinia)/CLOCKS_PER_SEC);
    std::cout << "Time for Newton: " << timea << std::endl;
    std::cout << "Number of DOF = " << cmeshwell->Solution().Rows() << std::endl;
    PostProcesswithVTK(an);
}

void TransientProblem(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh * > meshvec)
{
    TPZCompMesh * cmeshwell = CmeshMixed(gmesh);
    
    int wellmatid = 1;
    REAL day = 86400.0;
    REAL dt = 0.1*day*1.0e14;
    TPZMaterial * Matwell = cmeshwell->FindMaterial(wellmatid);
    TPZMonoPhaseWell *Well = dynamic_cast<TPZMonoPhaseWell *>(Matwell);
    Well->Setdt(dt);
    
    // Transferindo para a multifisica
    TPZBuildMultiphysicsMesh::AddElements(meshvec, cmeshwell);
    TPZBuildMultiphysicsMesh::AddConnects(meshvec, cmeshwell);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, cmeshwell);
    
#ifdef PZDEBUG
    std::ofstream dumpfile("ComputationaMeshMultiphysic.txt");
    cmeshwell->Print(dumpfile);
#endif
    
    
    // Setting up the analysis configuration
    int numofThreads = 4;
    TPZAnalysis * an = new TPZAnalysis(cmeshwell,true);
    TPZSkylineNSymStructMatrix skylnsym(cmeshwell);
    TPZStepSolver<STATE> step;
    skylnsym.SetNumThreads(numofThreads);
    step.SetDirect(ELU);
    an->SetStructuralMatrix(skylnsym);
    an->SetSolver(step);

    TimeForward(an, meshvec, cmeshwell);
}

void PrintLS(TPZAnalysis *an)
{
    TPZAutoPointer< TPZMatrix<STATE> > KGlobal;
    TPZFMatrix<STATE> FGlobal;
    KGlobal =   an->Solver().Matrix();
    FGlobal =   an->Rhs();
    
#ifdef PZDEBUG
#ifdef LOG4CXX
    if(logdata->isDebugEnabled())
    {
        std::stringstream sout;
        KGlobal->Print("KGlobal = ", sout,EMathematicaInput);
        FGlobal.Print("FGlobal = ", sout,EMathematicaInput);
        LOGPZ_DEBUG(logdata,sout.str())
    }
#endif
#endif
    
}

void PostProcesswithVTK(TPZAnalysis *an)
{
    const int dim = 1;
    int div = 0;
    TPZStack<std::string> scalnames, vecnames;
    std::string plotfile = "Well.vtk";
    
    scalnames.Push("Momentum");
    scalnames.Push("Pressure");
    scalnames.Push("Velocity");
    scalnames.Push("Density");
    an->DefineGraphMesh(dim, scalnames, vecnames, plotfile);
    an->PostProcess(div);
}

TPZCompMesh * CmeshMixed(TPZGeoMesh * gmesh)
{
    int dim = 1;
    int wellId = 1;
    int bottomId = 2;
    int topId = 3;
    int loadcases = 1;
    
    const int typeFlux = 1, typePressure = 0;
    TPZFMatrix<STATE> val1(3,2,0.), val2(3,1,0.);
    
    // Computational mesh
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    
    // Well Material
    TPZMonoPhaseWell * mat = new TPZMonoPhaseWell(wellId);
    mat->SetNumLoadCases(loadcases);
    cmesh->InsertMaterialObject(mat);
    
    // Rigth hand side function
    TPZDummyFunction<STATE> *dum = new TPZDummyFunction<STATE>(Ffunction, 5);
    TPZAutoPointer<TPZFunction<STATE> > forcef;
    dum->SetPolynomialOrder(0);
    forcef = dum;
    mat->SetForcingFunction(forcef);
    
    // Bc Bottom
    val2(0,0) = 100.0; //inlet Pressure
    TPZBndCond * bcBottom = mat->CreateBC(mat, bottomId, typeFlux, val1, val2);
    
    // Bc Top
    val2(0,0) = 10.0e6; // outlet Pressure
    TPZBndCond * bcTop = mat->CreateBC(mat, topId, typePressure, val1, val2);

    
    cmesh->InsertMaterialObject(bcBottom);
    cmesh->InsertMaterialObject(bcTop);
    
    cmesh->SetDimModel(dim);
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    cmesh->AutoBuild();
    
    return cmesh;
}

TPZCompMesh * CmeshMixedInitial(TPZGeoMesh * gmesh)
{
    int dim = 1;
    int wellId = 1;
    int bottomId = 2;
    int topId = 3;
    int loadcases = 1;
    
    const int typeFlux = 1, typePressure = 0;
    TPZFMatrix<STATE> val1(3,2,0.), val2(3,1,0.);
    
    // Computational mesh
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    
    // Well Material
    TPZMonoPhaseWell * mat = new TPZMonoPhaseWell(wellId);
    mat->SetNumLoadCases(loadcases);
    cmesh->InsertMaterialObject(mat);
    
    // Rigth hand side function
    TPZDummyFunction<STATE> *dum = new TPZDummyFunction<STATE>(Ffunction, 5);
    TPZAutoPointer<TPZFunction<STATE> > forcef;
    dum->SetPolynomialOrder(0);
    forcef = dum;
    mat->SetForcingFunction(forcef);
    
    // Bc Bottom
    val2(0,0) = 0.0; //inlet Pressure
    TPZBndCond * bcBottom = mat->CreateBC(mat, bottomId, typeFlux, val1, val2);
    
    // Bc Top
    val2(0,0) = 10.0e6; // outlet Pressure
    TPZBndCond * bcTop = mat->CreateBC(mat, topId, typePressure, val1, val2);
    
    
    cmesh->InsertMaterialObject(bcBottom);
    cmesh->InsertMaterialObject(bcTop);
    
    cmesh->SetDimModel(dim);
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    cmesh->AutoBuild();
    
    return cmesh;
}

TPZCompMesh * CmeshFlux(int qorder, TPZGeoMesh * gmesh)
{
    
    int dim = 1;
    int wellId = 1;
    int bottomId = 2;
    int topId = 3;
    int loadcases = 1;
    
    const int typeFlux = 0, typePressure = 1;
    TPZFMatrix<STATE> val1(3,2,0.), val2(3,1,0.);
    
    // Malha computacional
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    
    TPZMonoPhaseWell * mat = new TPZMonoPhaseWell(wellId);
    mat->SetNumLoadCases(loadcases);
    cmesh->InsertMaterialObject(mat);
    
    
    
    // Bc Bottom
    TPZBndCond * bcBottom = mat->CreateBC(mat, bottomId, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(bcBottom);
    
    // Bc Top
    TPZBndCond * bcTop = mat->CreateBC(mat, topId, typePressure, val1, val2);
    cmesh->InsertMaterialObject(bcTop);
    
    
    
    // Setando Hdiv
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(qorder);
    cmesh->SetAllCreateFunctionsContinuous();
    
    cmesh->AutoBuild();
    
    
#ifdef PZDEBUG
    std::ofstream out("cmeshFlux.txt");
    cmesh->Print(out);
#endif
    
    return cmesh;
}

TPZCompMesh * CmeshPressure(int porder, TPZGeoMesh * gmesh)
{
    
    int dim = 1;
    int wellId = 1;
    int bottomId = 2;
    int topId = 3;
    int loadcases = 1;
    
    const int typeFlux = 0, typePressure = 1;
    TPZFMatrix<STATE> val1(3,2,0.), val2(3,1,0.);
    
    // Computational mesh
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    
    TPZMonoPhaseWell * material = new TPZMonoPhaseWell(wellId);
    material->SetNumLoadCases(loadcases);
    cmesh->InsertMaterialObject(material);
    
    // Bc Bottom
    TPZMaterial * bcBottom = material->CreateBC(material, bottomId, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(bcBottom);
    
    // Bc Top
    TPZBndCond * bcTop = material->CreateBC(material, topId, typePressure, val1, val2);
    cmesh->InsertMaterialObject(bcTop);
    
    // Setando L2
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(porder);
    
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(false);
    
    cmesh->AutoBuild();
    
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
    
#ifdef PZDEBUG
    std::ofstream out("cmeshPress.txt");
    cmesh->Print(out);
#endif
    
    return cmesh;
    
}

void ParametricfunctionS(const TPZVec<REAL> &par, TPZVec<REAL> &X)
{
    X[0] = 0.0;
    X[1] = 0.0;
    X[2] = par[0];
}

void PrintGeoMesh(TPZGeoMesh * gmesh)
{
    
    
#ifdef PZDEBUG
    //  Print Geometrical Base Mesh
    std::ofstream argument("GeometicMesh.txt");
    gmesh->Print(argument);
    std::ofstream Dummyfile("GeometricMesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh,Dummyfile, true);
    
#endif
}

TPZGeoMesh * WellMesh(REAL s, REAL ds,int nelements)
{
    
    // Creating a 0D element to be extruded
    TPZGeoMesh * GeoMesh0D = new TPZGeoMesh;
    GeoMesh0D->NodeVec().Resize(1);
    TPZGeoNode Node;
    TPZVec<REAL> coors(3,0.0);
    Node.SetCoord(coors);
    Node.SetNodeId(0);
    GeoMesh0D->NodeVec()[0]=Node;
    
    TPZVec<int64_t> Topology(1,0);
    int elid=0;
    int matid=1;
    
    new TPZGeoElRefPattern < pzgeom::TPZGeoPoint >(elid,Topology,matid,*GeoMesh0D);
    GeoMesh0D->BuildConnectivity();
    GeoMesh0D->SetDimension(0);
    
    TPZHierarquicalGrid CreateGridFrom(GeoMesh0D);
    TPZAutoPointer<TPZFunction<REAL> > ParFuncX = new TPZDummyFunction<REAL>(ParametricfunctionS, 5);
    CreateGridFrom.SetParametricFunction(ParFuncX);
    CreateGridFrom.SetFrontBackMatId(2,3);
    
    // Computing Mesh extruded along the parametric curve -> ParametricfunctionS
    TPZGeoMesh * GeoMesh1D = CreateGridFrom.ComputeExtrusion(s, ds, nelements);
    
    return GeoMesh1D;
    
}

void AssembleLastState(TPZAnalysis *an, TPZCompMesh * cmesh){
    
    // Getting the material
    int wellmatid  = 1;
    bool state = false;
    TPZMaterial * Matwell = cmesh->FindMaterial(wellmatid);
    TPZMonoPhaseWell *Well = dynamic_cast<TPZMonoPhaseWell *>(Matwell);
    Well->SetNextStep(state);
    an->AssembleResidual();

    
}

void AssembleNextState(TPZAnalysis *an, TPZCompMesh * cmesh){

    // Getting the material
    int wellmatid  = 1;
    bool state = true;
    TPZMaterial * Matwell = cmesh->FindMaterial(wellmatid);
    TPZMonoPhaseWell *Well = dynamic_cast<TPZMonoPhaseWell *> (Matwell);
    Well->SetNextStep(state);
    an->Assemble();
    
}

void NewtonIterations(TPZAnalysis *an, TPZManVector<TPZCompMesh *> meshvector, TPZCompMesh * cmesh){
    
    TPZFMatrix<STATE> Residual(an->Rhs().Rows(),1,0.0);
    TPZFMatrix<STATE> ResidualAtn = Residual;
    TPZFMatrix<STATE> ResidualAtnplusOne = Residual;
    TPZFMatrix<STATE> alphaAtn = an->Solution();

    AssembleLastState(an,cmesh);
    ResidualAtn = an->Rhs();
    AssembleNextState(an,cmesh);
    ResidualAtnplusOne = an->Rhs();
    
    Residual = ResidualAtn + ResidualAtnplusOne;

    TPZFMatrix<STATE> X = alphaAtn;
    TPZFMatrix<STATE> DeltaX = alphaAtn;

    REAL error     =   1;
    REAL normdx    =   1;
    REAL toldx       =   1e-5;
    REAL tolrx       =   1e-5;
    
    int miterations  =   20;
    int iterations  =   0;
    int centinel    =   0;
    int fixed       =   0;
    
    while (error >= tolrx && iterations <= miterations) {

        an->Rhs() = Residual;
        an->Rhs() *= -1.0;

        an->Solve();
        
        DeltaX = an->Solution();
        normdx = Norm(DeltaX);
        X += DeltaX;

        cmesh->LoadSolution(X);
        TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvector, cmesh);

        if (((fixed+1) * (centinel) == iterations)) {
            an->Assemble();
            centinel++;
        }
        else{
            
            an->AssembleResidual();
        }
        
        ResidualAtnplusOne = an->Rhs();
        Residual = ResidualAtn + ResidualAtnplusOne;
        error = Norm(Residual);
        iterations++;
        
#ifdef PZDEBUG
    #ifdef LOG4CXX
            if(logdata->isDebugEnabled())
            {
                std::stringstream sout;
                ResidualAtn.Print("ResidualAtn = ", sout,EMathematicaInput);
                ResidualAtnplusOne.Print("ResidualAtnplusOne = ", sout,EMathematicaInput);
                DeltaX.Print("DeltaX = ", sout,EMathematicaInput);
                X.Print("X = ", sout,EMathematicaInput);
                LOGPZ_DEBUG(logdata,sout.str())
            }
    #endif
#endif

        if(error < tolrx || normdx < toldx)
        {
            std::cout << "Converged with iterations:  " << iterations << std::endl;
            std::cout << "error norm: " << error << std::endl;
            std::cout << "error of dx: " << normdx << std::endl;
            break;
        }

        if (iterations == miterations) {
            std::cout << "Out max iterations " << iterations << std::endl;
            std::cout << "error norm " << error << std::endl;
            break;
        }
        
    }
    
    
}

void TimeForward(TPZAnalysis *an, TPZManVector<TPZCompMesh *> meshvector, TPZCompMesh * cmesh){
    
    REAL day = 86400.0;
    REAL maxtime = 10.0*day;
    REAL dt = 0.1*day;
    
    int nsteps = maxtime / dt;
    REAL tk = 0;

    for (int istep = 1 ; istep <=  nsteps; istep++) {
        tk = istep*dt;

        std::cout << "Time Value (days): " << tk/86400.0 << std::endl;
        std::cout<<  "Time step: " << istep << std::endl;

        const clock_t tinia = clock();
        NewtonIterations(an,meshvector,cmesh);
        const clock_t tenda = clock();
        const REAL timea = REAL(REAL(tenda - tinia)/CLOCKS_PER_SEC);
        std::cout << "Time for Newton: " << timea << std::endl;
        std::cout << "Number of DOF = " << cmesh->Solution().Rows() << std::endl;
        PostProcesswithVTK(an);
        
        
        
    }
    
    
}

void Ffunction(const TPZVec<REAL> &pt, TPZVec<STATE> &ff)
{    
//    if (fabs(z-100.0) <= 1.0 ) {
//        ff[0] = 10000.0;
//    }
//    
//    ff[0] = 100.0;
}

