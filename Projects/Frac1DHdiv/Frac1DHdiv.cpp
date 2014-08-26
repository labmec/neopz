#include "pzvec.h"
#include "pzadmchunk.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "TPZInterfaceEl.h"
#include "pzgnode.h"
#include "pzmaterial.h"
#include "pzfmatrix.h"
#include "pzerror.h"
#include "pzgeoel.h"
#include "pzmatrix.h"
#include "TPZReadGIDGrid.h"
#include "TPZVTKGeoMesh.h"
#include "time.h"
#include "pzconvectionproblem.h"
#include "TPZMatfrac1dhdiv.h"
#include "pzl2projection.h"
#include "TPZSkylineNSymStructMatrix.h"

#include "pzmultiphysicscompel.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzmultiphysicselement.h"
#include "TPZMultiphysicsInterfaceEl.h"

#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "pzstepsolver.h"
#include "pzintel.h"
#include "pzbndcond.h"
#include "pzlog.h"

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "pzgeopoint.h"
#include "TPZGeoLinear.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"
#include "pzgeoelside.h"

#include "pzequationfilter.h"
#include "pzgradientreconstruction.h"
#include "pzl2projection.h"
#include "pzbfilestream.h"

#include "pzskylstrmatrix.h"
#include "pzfstrmatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzpoisson3d.h"
#include <time.h>
#include <stdio.h>

// Using Log4cXX as logging tool
#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.multiphase"));
#endif

#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.material.multiphase.data"));
#endif
//
// End Using Log4cXX as logging tool

TPZGeoMesh *LineGeometry();

TPZCompMesh *ComputationalMeshBulkflux(TPZGeoMesh *gmesh, int pOrder);
TPZCompMesh *ComputationalMeshPseudopressure(TPZGeoMesh *gmesh, int pOrder);
TPZCompMesh *ComputationalMeshMixed(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec);
TPZCompMesh *L2ProjectionQ(TPZGeoMesh *gmesh, int pOrder, TPZVec<STATE> &solini);
TPZCompMesh *L2ProjectionP(TPZGeoMesh *gmesh, int pOrder, TPZVec<STATE> &solini);
void RotateGeomesh(TPZGeoMesh *gmesh, REAL CounterClockwiseAngle);

void SolveSyst(TPZAnalysis &an, TPZCompMesh *cmesh);
void InitialFlux(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void InitialPressure(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void PosProcessBulkflux(TPZAnalysis &an, std::string plotfile);
void PosProcessL2(TPZAnalysis &an, std::string plotfile);
void PosProcessMultphysics(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile);
void UniformRefinement(TPZGeoMesh *gMesh, int nh);
void RefinUniformElemComp(TPZCompMesh  *cMesh, int ndiv);
void SolveSystemTransient(REAL deltaT,REAL maxTime, TPZAnalysis *NonLinearAn, TPZAnalysis *NonLinearAnTan, TPZVec<TPZCompMesh *> meshvec,TPZCompMesh* mphysics);
void BulkFlux(const TPZVec<REAL> &pt, TPZVec<STATE> &Val);

bool ftriang = false;
REAL angle = 0.0*M_PI/4.0;

int main()
{
	std::string dirname = PZSOURCEDIR;
#ifdef LOG4CXX
    std::string FileName = dirname;
    //   FileName = dirname + "/Projects/OilWaterSystem/";
    //   FileName += "OilWaterLog4cxx.cfg";
    //   InitializePZLOG(FileName);
    //   InitializePZLOG();
#endif
    
    //
    //#ifdef LOG4CXX
    //  InitializePZLOG();
    //#endif
    
    //      gRefDBase.InitializeAllUniformRefPatterns();
    
    //  Reading mesh
    std::string GridFileName;
    GridFileName = dirname + "/Projects/OilWaterSystem/";
    //    GridFileName += "OilWaterSystemUnit.dump";
    //    GridFileName += "Labyrinth.dump";
    //    GridFileName += "BaseGeometryMazeOne.dump";
    GridFileName += "LinearMesh.dump";
    
    //  GridFileName = "Labyrinth.dump";
    //  GridFileName = "OilWaterSystemUnitTwo.dump";
    //  GridFileName = "OilWaterSystemUnitOneHRef.dump";
    //  GridFileName = "OilWaterSystemUnitTwoHRef.dump";
    
//    TPZReadGIDGrid GeometryInfo;
//    GeometryInfo.SetfDimensionlessL(1.0);
//    TPZGeoMesh * gmesh = GeometryInfo.GeometricGIDMesh(GridFileName);
    TPZGeoMesh * gmesh = LineGeometry();
    RotateGeomesh(gmesh, angle);
    {
        //  Print Geometrical Base Mesh
        std::ofstream argument("GeometicMesh.txt");
        gmesh->Print(argument);
        std::ofstream Dummyfile("GeometricMesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh,Dummyfile, true);
    }
    
    //  std::vector<REAL> dd(2,0);
    
    
    int Href = 0;
    int div = 0;
    int POrderBulkFlux = 1;
    int POrderPseudopressure = 1;
    
    UniformRefinement(gmesh, Href);
    
    {
        //  Print Geometrical Base Mesh
        std::ofstream argument("RefGeometicMesh.txt");
        gmesh->Print(argument);
        std::ofstream Dummyfile("RefGeometricMesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh,Dummyfile, true);
    }
    
    
    // Computational meshes
    
    //  First computational mesh
    TPZCompMesh * CMeshBulkflux = ComputationalMeshBulkflux(gmesh, POrderBulkFlux);
    //  Print Second computational mesh
    std::ofstream ArgumentBulkflux("ComputationalMeshForBulkflux.txt");
    CMeshBulkflux->Print(ArgumentBulkflux);
    
    //  Second computational mesh
    TPZCompMesh * CMeshPseudopressure = ComputationalMeshPseudopressure(gmesh, POrderPseudopressure);
    //  Print First computational mesh
    std::ofstream ArgumentPseudopressure("ComputationalMeshForPseudopressure.txt");
    CMeshPseudopressure->Print(ArgumentPseudopressure);
    
    TPZAnalysis Anbulkflux(CMeshBulkflux);
    std::string outputfile1;
    outputfile1 = "SolutionBulkflux";
    std::stringstream outputfiletemp1;
    outputfiletemp1 << outputfile1 << ".vtk";
    std::string plotfilebuklflux = outputfiletemp1.str();
    TPZFMatrix<STATE> InitialQSolution = Anbulkflux.Solution();
    int rwosQ= InitialQSolution.Rows();
    
    Anbulkflux.LoadSolution(InitialQSolution);
    int num= InitialQSolution.Rows();
    
    
    TPZVec<STATE> soliniQ(num,0.0);
    TPZCompMesh  * cmeshQL2 = L2ProjectionQ(gmesh, POrderBulkFlux, soliniQ);
    
    TPZAnalysis anQL2(cmeshQL2);
    //    SolveSyst(anQL2, cmeshQL2);
    //    Anbulkflux.LoadSolution(InitialQSolution);
    PosProcessBulkflux(Anbulkflux,plotfilebuklflux);
    
    
    
    TPZAnalysis AnPressure(CMeshPseudopressure);
    
    std::string outputfile2;
    outputfile2 = "SolutionPressure";
    std::stringstream outputfiletemp2;
    outputfiletemp2 << outputfile2 << ".vtk";
    std::string plotfilePressure = outputfiletemp2.str();
    TPZFMatrix<STATE> InitialPSolution = AnPressure.Solution();
    //  for (int i=0; i < InitialPSolution.Rows(); i++) {
    //      InitialPSolution(i)=0.0;
    //  }
    
    
    TPZVec<STATE> solini(InitialPSolution.Rows(),0.0);
    TPZCompMesh  * cmeshL2 = L2ProjectionP(gmesh, POrderPseudopressure, solini);
    
    TPZAnalysis anL2(cmeshL2);
    SolveSyst(anL2, cmeshL2);
    
    AnPressure.LoadSolution(anL2.Solution());
    PosProcessL2(AnPressure,plotfilePressure);
    
    
    //  Multiphysics Mesh
    TPZVec<TPZCompMesh *> meshvec(2);//4);
    meshvec[0] = CMeshBulkflux;
    meshvec[1] = CMeshPseudopressure;
    
    
    TPZCompMesh * MultiphysicsMesh = ComputationalMeshMixed(gmesh,meshvec);
    std::ofstream ArgumentMultiphysic("MultiphysicsMesh.txt");
    MultiphysicsMesh->Print(ArgumentMultiphysic);
    
    
    TPZAnalysis *MultiphysicsAn = new TPZAnalysis(MultiphysicsMesh);
    TPZAnalysis *MultiphysicsAnTan = new TPZAnalysis(MultiphysicsMesh);
#ifdef DEBUG
    int Nthreads = 1;
#else
    int Nthreads = 4;
#endif
    
    TPZSkylineNSymStructMatrix matsk(MultiphysicsMesh);
    //TPZSkylineNSymStructMatrix matskTan(MultiphysicsMesh);
    
    //TPZFrontStructMatrix <TPZFrontNonSym<STATE> > *matsk = new TPZFrontStructMatrix <TPZFrontNonSym<STATE> > (MultiphysicsMesh);
    
    //matsk->SetQuiet(0);
    MultiphysicsAn->SetStructuralMatrix(matsk);
    MultiphysicsAn->StructMatrix()->SetNumThreads(Nthreads);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELU);
    MultiphysicsAn->SetSolver(step);
    
    
    std::string outputfile;
    outputfile = "TransientSolutionini";
    std::stringstream outputfiletemp;
    outputfiletemp << outputfile << ".vtk";
    std::string plotfile = outputfiletemp.str();
    PosProcessMultphysics(meshvec,MultiphysicsMesh,*MultiphysicsAn,plotfile);
    
    //     Tima control parameters
    
    REAL hour = 60.0*60.0;
    REAL day = 24.0*hour;
    REAL year = 365.0*day;
    
    REAL deltaT = 0.0025;
    REAL maxTime = 0.5;
    SolveSystemTransient(deltaT, maxTime, MultiphysicsAn, MultiphysicsAnTan, meshvec, MultiphysicsMesh);
    return 0;
    
}


TPZCompMesh *ComputationalMeshBulkflux(TPZGeoMesh *gmesh, int pOrder)
{
    /// criar materiais
    int dim = 1;
    int matId1 = 1;
    
    TPZMatPoisson3d *material1;
    material1 = new TPZMatPoisson3d(matId1,dim);
    TPZMaterial * mat1(material1);
    material1->NStateVariables();
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
    cmesh->InsertMaterialObject(mat1);
    
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    TPZMaterial * BCond2 = material1->CreateBC(mat1,2,0, val1, val2);
    TPZMaterial * BCond3 = material1->CreateBC(mat1,3,0, val1, val2);
    
    cmesh->SetAllCreateFunctionsHDiv();// Hdiv approximation space
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    
    
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    return cmesh;
}


TPZCompMesh *ComputationalMeshPseudopressure(TPZGeoMesh *gmesh, int pOrder)
{
    /// criar materiais
    int dim = 1;
    int matId1 = 1;
    
    TPZMatPoisson3d *material1;
    material1 = new TPZMatPoisson3d(matId1,dim);
    TPZMaterial * mat1(material1);
    material1->NStateVariables();

    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
    cmesh->InsertMaterialObject(mat1);
    
    
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    TPZMaterial * BCond2 = material1->CreateBC(mat1,2,0, val1, val2);
    TPZMaterial * BCond3 = material1->CreateBC(mat1,3,0, val1, val2);
    
    cmesh->SetAllCreateFunctionsDiscontinuous(); // L2 approximation space
    
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    
    
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
    cmesh->AutoBuild();
    
    
    ///inserir connect da pressao
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        //newnod.SetPressure(true);
        newnod.SetLagrangeMultiplier(1);
    }
    
    ///set order total da shape
    int nel = cmesh->NElements();
    for(int i=0; i<nel; i++){
        TPZCompEl *cel = cmesh->ElementVec()[i];
        TPZCompElDisc *celdisc = dynamic_cast<TPZCompElDisc *>(cel);
        celdisc->SetConstC(1.);
        celdisc->SetCenterPoint(0, 0.);
        celdisc->SetCenterPoint(1, 0.);
        celdisc->SetCenterPoint(2, 0.);
        celdisc->SetTrueUseQsiEta();
        if(celdisc && celdisc->Reference()->Dimension() == cmesh->Dimension())
        {
            if(ftriang==true || celdisc->Reference()->Type()==ETriangle) celdisc->SetTotalOrderShape();
            else celdisc->SetTensorialShape();
        }
    }
    
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    
    return cmesh;
}

TPZCompMesh *ComputationalMeshMixed(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec){
    
    
    //Creating computational mesh for multiphysic elements
    gmesh->ResetReference();
    TPZCompMesh *mphysics = new TPZCompMesh(gmesh);;
    
    bool Dimensionless = true;
    
    int dim =1;
    int matId1 = 1;
    // Setting data
    
    REAL viscosity = 1.;
    TPZMatfrac1dhdiv *material1 = new TPZMatfrac1dhdiv(matId1,dim,viscosity);
    
    REAL deltaT = 0.1;
    REAL maxTime = 0.1;
    REAL MPa = 1.0e+6;
    
    if (Dimensionless)
    {
        REAL Rhoref,Etaref,Lref,Kref,Pref,Qref;
        Rhoref = 1000.0;
        Etaref = 1.0e-3;
        Lref = 100.0;
        Kref = 1.0e-13;
        Pref = 20.0*MPa;
        Qref = (Rhoref*(Kref/Etaref))*(Pref/Lref);
        
    }
    
    std::string GridFileName, dirname = PZSOURCEDIR;
    GridFileName = dirname + "/Projects/OilWaterSystem/";
    GridFileName += "LabyrinthKvalues.txt";
    
    material1->SetTimeStep(1.0);
    material1->SetTime(0.0);
    material1->SetTScheme(1.0);
    
    TPZMaterial *mat1(material1);
    mphysics->InsertMaterialObject(mat1);
    mphysics->SetDimModel(dim);
    
    
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);
    
    val2(0,0)=1.0*0.0020*cos(angle);// qx
    val2(1,0)=0.0*20.0*MPa;// P
    TPZMaterial * BCondL = material1->CreateBC(mat1,2,0, val1, val2);
    
    
    val2(0,0)=0.000;// qx
    val2(1,0)=1.0*18.0*MPa;// P
    TPZMaterial * BCondR = material1->CreateBC(mat1,3,1, val1, val2);
    
    
    mphysics->SetAllCreateFunctionsMultiphysicElem();
    mphysics->InsertMaterialObject(BCondL);
    mphysics->InsertMaterialObject(BCondR);
    
    
    mphysics->AutoBuild();
    mphysics->AdjustBoundaryElements();
    mphysics->CleanUpUnconnectedNodes();
    
    // Creating multiphysic elements into mphysics computational mesh
    TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
    TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
    
    mphysics->Reference()->ResetReference();
    mphysics->LoadReferences();
    
    // Creation of interface elements
    int nel = mphysics->ElementVec().NElements();
    for(int el = 0; el < nel; el++)
    {
        TPZCompEl * compEl = mphysics->ElementVec()[el];
        if(!compEl) continue;
        int index = compEl ->Index();
        if(compEl->Dimension() == mphysics->Dimension())
        {
            TPZMultiphysicsElement * InterpEl = dynamic_cast<TPZMultiphysicsElement *>(mphysics->ElementVec()[index]);
            if(!InterpEl) continue;
            InterpEl->CreateInterfaces();
        }
    }
    
    return mphysics;
}

void PosProcessBulkflux(TPZAnalysis &an, std::string plotfile){
    TPZManVector<std::string,10> scalnames(0), vecnames(1);
    vecnames[0]= "FluxL2";
    
    const int dim = 2;
    int div = 0;
    an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
    an.PostProcess(div,dim);
    std::ofstream out("malhaflux.txt");
    an.Print("nothing",out);
}

void PosProcessL2(TPZAnalysis &an, std::string plotfile){
    TPZManVector<std::string,10> scalnames(1), vecnames(0);
    scalnames[0]= "Solution";
    
    const int dim = 2;
    int div = 0;
    an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
    an.PostProcess(div,dim);
    std::ofstream out("malhaflux.txt");
    an.Print("nothing",out);
}

void PosProcessMultphysics(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile){
    
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
    TPZManVector<std::string,10> scalnames(1), vecnames(1);
    
    scalnames[0] = "Pressure";
    vecnames[0] = "BulkVelocity";
    
    const int dim = 2;
    int div =1;
    an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
    an.PostProcess(div,dim);
    std::ofstream out("malha.txt");
    an.Print("nothing",out);
    
}

void UniformRefinement(TPZGeoMesh *gMesh, int nh)
{
    for ( int ref = 0; ref < nh; ref++ ){
        TPZVec<TPZGeoEl *> filhos;
        long n = gMesh->NElements();
        for ( long i = 0; i < n; i++ ){
            TPZGeoEl * gel = gMesh->ElementVec() [i];
            if (gel->Dimension() == 2 || gel->Dimension() == 1) gel->Divide (filhos);
        }//for i
    }//ref
}

void RefinUniformElemComp(TPZCompMesh  *cMesh, int ndiv)
{
    
    TPZVec<long > subindex;
    for (long iref = 0; iref < ndiv; iref++) {
        TPZAdmChunkVector<TPZCompEl *> elvec = cMesh->ElementVec();
        long nel = elvec.NElements();
        for(long el=0; el < nel; el++){
            TPZCompEl * compEl = elvec[el];
            if(!compEl) continue;
            long ind = compEl->Index();
            compEl->Divide(ind, subindex, 0);
        }
    }
    
}

void SolveSystemTransient(REAL deltaT,REAL maxTime, TPZAnalysis *NonLinearAn, TPZAnalysis *NonLinearAnTan, TPZVec<TPZCompMesh *> meshvec,TPZCompMesh* mphysics){
    
    TPZFMatrix<STATE> SolutiontoLoad;
    TPZFMatrix<STATE> SolutiontoSave = meshvec[2]->Solution();
    //    RotateGeomesh(mphysics->Reference(),angle);
    //  {
    //      TPZBFileStream load;
    //      load.OpenRead("MultiphaseSaturationSol.bin");
    //      SolutiontoLoad.Read(load,0);
    //      meshvec[2]->LoadSolution(SolutiontoLoad);
    //      TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);
    //  }
    
    std::string OutPutFile = "TransientSolution";
    TPZMaterial *mat1 = mphysics->FindMaterial(1);
    //    TPZMaterial *mat2 = mphysics->FindMaterial(2);
    
    TPZMatfrac1dhdiv * material1 = dynamic_cast<TPZMatfrac1dhdiv *>(mat1);

    material1->SetTimeStep(deltaT);
    material1->SetTime(0.0);
    material1->SetTScheme(1.0);
    
    //  Starting Newton Iterations
    TPZFMatrix<STATE> DeltaX = mphysics->Solution();
    TPZFMatrix<STATE> Uatn = mphysics->Solution();
    TPZFMatrix<STATE> Uatk = mphysics->Solution();
    
    
    REAL TimeValue = 0.0;
    REAL Tolerance = 1.0e-8;
    int cent = 0;
    int MaxIterations = 50;
    TimeValue = cent*deltaT;
    REAL NormValue =1.0;
    bool StopCriteria = false;
    TPZFMatrix<STATE> RhsAtn, RhsAtnPlusOne, Residual;
    TPZFMatrix<STATE> RhsAtnT, RhsAtnPlusOneT, ResidualT;
    
    std::string outputfile;
    outputfile = OutPutFile;
    std::stringstream outputfiletemp;
    outputfiletemp << outputfile << ".vtk";
    std::string plotfile = outputfiletemp.str();
    PosProcessMultphysics(meshvec,mphysics,*NonLinearAn,plotfile);
    
    
	TPZManVector<STATE> PrintStep(5);
	int control = 0;
	PrintStep[0]=0.1;
	PrintStep[1]=0.2;
	PrintStep[2]=0.3;
	PrintStep[3]=0.4;
	PrintStep[4]=0.5;
    
    
    material1->SetCurrentState();
    NonLinearAn->Assemble();
    RhsAtnPlusOne = NonLinearAn->Rhs();
    
    
    // Starting
    material1->SetCurrentState();
    NonLinearAn->Assemble();
    RhsAtnPlusOne = NonLinearAn->Rhs();
    
    std::cout << " Starting the time computations. " << std::endl;
    std::cout << " Number of DOF: " << mphysics->Solution().Rows() <<  std::endl;
    while (TimeValue < maxTime)
    {
        
        
        material1->SetLastState();
        NonLinearAn->AssembleResidual();
        RhsAtn = NonLinearAn->Rhs();
        
        //         material1->SetCurrentState();
        //         NonLinearAn->Assemble();
        //         RhsAtnPlusOne = NonLinearAn->Rhs();
        Residual= RhsAtn + RhsAtnPlusOne;
        NormValue = Norm(Residual);
        
        //      material1->SetLastState();
        //      NonLinearAnTan->Assemble();
        //      RhsAtnT = NonLinearAnTan->Rhs();
        //
        //      material1->SetCurrentState();
        //      NonLinearAnTan->Assemble();
        //      RhsAtnPlusOneT = NonLinearAnTan->Rhs();
        //      ResidualT= RhsAtnT + RhsAtnPlusOneT;
        //      NormValueT = Norm(ResidualT);
        
#ifdef LOG4CXX
        if(logdata->isDebugEnabled())
        {
            std::stringstream sout;
            Residual.Print("Residual = ",sout,EMathematicaInput);
            LOGPZ_DEBUG(logdata,sout.str());
        }
#endif
        
        TPZAutoPointer< TPZMatrix<REAL> > matK;
        TPZFMatrix<STATE> fvec;
        matK=NonLinearAn->Solver().Matrix();
        fvec = NonLinearAn->Rhs();
        /*#ifdef LOG4CXX
         if(logdata->isDebugEnabled())
         {
         std::stringstream sout;
         matK->Print("matK = ", sout,EMathematicaInput);
         fvec.Print("fvec = ", sout,EMathematicaInput);
         LOGPZ_DEBUG(logdata,sout.str())
         }
         #endif*/
        
        int iterations= 0;
        while (NormValue > Tolerance)
        {
            
            Residual*=-1.0;
            NonLinearAn->Rhs()=Residual;
            const clock_t tini = clock();
            NonLinearAn->Solve();
            const clock_t tend = clock();
            const REAL time = REAL(REAL(tend - tini)/CLOCKS_PER_SEC);
            std::cout << "Time for solving: " << time << std::endl;
            DeltaX = NonLinearAn->Solution();
            Uatk = (Uatn + DeltaX);
            
            
            mphysics->LoadSolution(Uatn + DeltaX);
            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
            
#ifdef LOG4CXX
            if(logdata->isDebugEnabled())
            {
                std::stringstream sout;
                sout.precision(20);
                Residual.Print("Residual = ",sout,EMathematicaInput);
                Uatk.Print("Uatk = ",sout,EMathematicaInput);
                LOGPZ_DEBUG(logdata,sout.str());
            }
#endif
            
            //#ifdef LOG4CXX
            //          if(logdata->isDebugEnabled())
            //          {
            //              std::stringstream sout;
            //              Uatk.Print("Uatk = ",sout,EMathematicaInput);
            //              LOGPZ_DEBUG(logdata,sout.str());
            //          }
            //#endif
            
            
    
            
            
            //          material1->SetCurrentState();
            //          NonLinearAnTan->Assemble();
            //          RhsAtnPlusOneT = NonLinearAnTan->Rhs();
            //          ResidualT= RhsAtnT + RhsAtnPlusOneT;
            //          NormValueT = Norm(ResidualT);
            
            material1->SetCurrentState();
            material1->SetTime(cent*deltaT);
            const clock_t tini2 = clock();
            NonLinearAn->Assemble();
            RhsAtnPlusOne = NonLinearAn->Rhs();
            Residual= RhsAtn + RhsAtnPlusOne;
            NormValue = Norm(Residual);
            const clock_t tend2 = clock();
            const REAL time2 = REAL(REAL(tend2 - tini2)/CLOCKS_PER_SEC);
            std::cout << "Time for Assemble and computing norm: " << time2 << std::endl;
            
            
#ifdef LOG4CXX
            if(logdata->isDebugEnabled())
            {
                std::stringstream sout;
                sout.precision(15);
                Uatk.Print(sout);
                Residual.Print("Res = ",sout,EMathematicaInput);
                LOGPZ_DEBUG(logdata,sout.str());
            }
#endif
            
            
            //#ifdef LOG4CXX
            //          if(logdata->isDebugEnabled())
            //          {
            //              std::stringstream sout;
            //              sout.precision(20);
            //              Residual.Print(sout);
            //              NonLinearAn->Solution().Print(sout);
            //              LOGPZ_DEBUG(logdata,sout.str());
            //          }
            //#endif
            
            //          NonLinearAn->StructMatrix()->EquationFilter().Reset();
            //          NonLinearAn->StructMatrix()->EquationFilter().SetActiveEquations(NoGradients);
            //
            //          numofequactive = NonLinearAn->StructMatrix()->EquationFilter().NActiveEquations();
            
            iterations++;
            std::cout << " Newton's Iteration = : " << iterations  << "     L2 norm = : " << NormValue <<  std::endl;
            if (iterations == MaxIterations)
            {
                StopCriteria = true;
                std::cout << " Time Step number = : " << iterations  << "\n Exceed max iterations numbers = : " << MaxIterations <<  std::endl;
                break;
            }
            
            
            Uatn = Uatk;
            
        }
        
        
        //      TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
        //      SolutiontoSave = mphysics->Solution();
        //      TPZBFileStream save;
        //        save.OpenWrite("GRSolution.bin");
        //      SolutiontoSave.Write(save,0);
        
		if (fabs(PrintStep[control] - TimeValue) < 1.0e-8 || PrintStep.size()-1 == control)
		{
            
            
			const clock_t tini3 = clock();
			
			outputfile = OutPutFile;
			std::stringstream outputfiletemp;
			outputfiletemp << outputfile << ".vtk";
			std::string plotfile = outputfiletemp.str();
			PosProcessMultphysics(meshvec,mphysics,*NonLinearAn,plotfile);
			
			const clock_t tend3 = clock();
			const REAL time3 = REAL(REAL(tend3 - tini3)/CLOCKS_PER_SEC);
			std::cout << "Time for printing: " << time3 << std::endl;
			std::cout << "Control: " << control << std::endl;
			control++;
            
		}
        
        if (StopCriteria) {
            std::cout << " Newton's Iteration = : " << iterations  << "     L2 norm = : " << NormValue <<  std::endl;
            break;
        }
        
        cent++;
        TimeValue = cent*deltaT;
        
        std::cout << " Time Step :  " << cent  << "  Time :  " << TimeValue <<  std::endl;
        
        /*        if(cent==15)
         {
         UsingGradient=true;
         std::cout << " Now Using Gradient Reconstruction. " << std::endl;
         }  */
        
    }
    
    //  CheckConvergence(RhsAtn,NonLinearAn, meshvec, mphysics);
    //  CheckElConvergence(RhsAtn,NonLinearAn, meshvec, mphysics);
    
    
}


// Setting up initial conditions

void SolveSyst(TPZAnalysis &an, TPZCompMesh *Cmesh)
{
    
    TPZSkylineStructMatrix full(Cmesh);
    an.SetStructuralMatrix(full);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    //step.SetDirect(ELU);
    an.SetSolver(step);
    an.Run();
    
    //  //Saida de Dados: solucao e  grafico no VT
    //  ofstream file("Solutout");
    //  an.Solution().Print("solution", file);    //Solution visualization on Paraview (VTK)
}

TPZCompMesh *L2ProjectionP(TPZGeoMesh *gmesh, int pOrder, TPZVec<STATE> &solini)
{
    /// criar materiais
    int dim = 2;
    TPZL2Projection *material;
    material = new TPZL2Projection(1, dim, 1, solini, pOrder);
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    TPZMaterial * mat(material);
    cmesh->InsertMaterialObject(mat);
    TPZAutoPointer<TPZFunction<STATE> > forcef = new TPZDummyFunction<STATE>(InitialPressure);
    material->SetForcingFunction(forcef);
    cmesh->SetAllCreateFunctionsDiscontinuous();
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
    cmesh->AutoBuild();
    
    ///set order total da shape
    int nel = cmesh->NElements();
    for(int i=0; i<nel; i++){
        TPZCompEl *cel = cmesh->ElementVec()[i];
        TPZCompElDisc *celdisc = dynamic_cast<TPZCompElDisc *>(cel);
        if(celdisc && celdisc->Reference()->Dimension() == cmesh->Dimension())
        {
            if(ftriang==true || celdisc->Reference()->Type()==ETriangle) celdisc->SetTotalOrderShape();
            else celdisc->SetTensorialShape();
        }
    }
    
    return cmesh;
    
}

TPZCompMesh *L2ProjectionQ(TPZGeoMesh *gmesh, int pOrder, TPZVec<STATE> &solini)
{
    /// criar materiais
    int dim = 2;
    TPZL2Projection *material;
    material = new TPZL2Projection(1, dim, 1, solini, pOrder);
    
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(dim);
    TPZMaterial * mat(material);
    cmesh->InsertMaterialObject(mat);
    TPZAutoPointer<TPZFunction<STATE> > forcef = new TPZDummyFunction<STATE>(InitialFlux);
    material->SetForcingFunction(forcef);
    cmesh->SetAllCreateFunctionsDiscontinuous();
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
    cmesh->AutoBuild();
    
    return cmesh;
    
}


// It requires modfify L2 number os state variables
void InitialFlux(const TPZVec<REAL> &pt, TPZVec<STATE> &disp)
{
    REAL x = pt[0];
    REAL y = pt[1];
    disp[0] = 0.0;
    //    disp[1] = 0.0;
    
}

void InitialPressure(const TPZVec<REAL> &pt, TPZVec<STATE> &disp)
{
    REAL x = pt[0];
    REAL y = pt[1];
    disp[0] = 0.0*(1.0 - 0.1 * x);
    
}












void RotateGeomesh(TPZGeoMesh *gmesh, REAL CounterClockwiseAngle)
{
    REAL theta = CounterClockwiseAngle;
    // It represents a 3D rotation around the z axis.
    TPZFMatrix<STATE> RotationMatrix(3,3,0.0);
    RotationMatrix(0,0) =   +cos(theta);
    RotationMatrix(0,1) =   -sin(theta);
    RotationMatrix(1,0) =   +sin(theta);
    RotationMatrix(1,1) =   +cos(theta);
    RotationMatrix(2,2) = 1.0;
    TPZVec<STATE> iCoords(3,0.0);
    TPZVec<STATE> iCoordsRotated(3,0.0);
    
    RotationMatrix.Print("Rotation = ");
    
    int NumberofGeoNodes = gmesh->NNodes();
    for (int inode = 0; inode < NumberofGeoNodes; inode++)
    {
        TPZGeoNode GeoNode = gmesh->NodeVec()[inode];
        GeoNode.GetCoordinates(iCoords);
        // Apply rotation
        iCoordsRotated[0] = RotationMatrix(0,0)*iCoords[0]+RotationMatrix(0,1)*iCoords[1]+RotationMatrix(0,2)*iCoords[2];
        iCoordsRotated[1] = RotationMatrix(1,0)*iCoords[0]+RotationMatrix(1,1)*iCoords[1]+RotationMatrix(1,2)*iCoords[2];
        iCoordsRotated[2] = RotationMatrix(2,0)*iCoords[0]+RotationMatrix(2,1)*iCoords[1]+RotationMatrix(2,2)*iCoords[2];
        GeoNode.SetCoord(iCoordsRotated);
        gmesh->NodeVec()[inode] = GeoNode;
    }
    
}

TPZGeoMesh *LineGeometry(){
    
    int Nodes = 2;
    int matid= 1;
    int bcl = 2;
    int bcr = 3;
	
	TPZGeoMesh * gmesh = new TPZGeoMesh;
	gmesh->SetMaxNodeId(Nodes);
	gmesh->NodeVec().Resize(Nodes);
	TPZVec<TPZGeoNode> Node(Nodes);
	
	TPZVec <long> TopolLine(2);
    TPZVec <long> TopolPoint(1);
	
	//indice dos nos
	long id = 0;

    Node[id].SetNodeId(id);
    Node[id].SetCoord(0 ,0.0 );//coord X
    Node[id].SetCoord(1 ,0.0 );//coord Y
    gmesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    Node[id].SetCoord(0 ,1.0 );//coord X
    Node[id].SetCoord(1 ,0.0 );//coord Y
    gmesh->NodeVec()[id] = Node[id];
    id++;
    
	
	
	//indice dos elementos
	id = 0;
    TopolLine[0] = 0;
    TopolLine[1] = 1;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,matid,*gmesh);
    id++;
    
    TopolPoint[0] = 0;
    new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (id,TopolPoint,bcl,*gmesh);
    id++;
    
    TopolPoint[0] = 1;
    new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (id,TopolPoint,bcr,*gmesh);
    id++;
    
	gmesh->BuildConnectivity();
    
    
	return gmesh;
}

