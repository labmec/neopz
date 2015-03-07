/*
 *  pznondarcyanalysis.cpp
 *  PZ
 *
 *  Created by Omar Duran Triana on 5/21/14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */

#include "TPZDarcyAnalysis.h"

#include <pzlog.h>
#include "TPZReadGIDGrid.h"
#include "TPZVTKGeoMesh.h"
#include "TPZAxiSymmetricDarcyFlow.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "math.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.DarcyFlow"));
#endif

TPZDarcyAnalysis::TPZDarcyAnalysis(TPZAutoPointer<SimulationData> DataSimulation,TPZVec<TPZAutoPointer<ReservoirData> > Layers)
{

    fgmesh=NULL;
    fmeshvec.Resize(2);
    fcmeshMixed=NULL;
    fResidualLastState.Resize(1, 1);
    fResidualLastState.Zero();
    fSimulationData = DataSimulation;
    fLayers = Layers;

}


TPZDarcyAnalysis::~TPZDarcyAnalysis()
{
	
}

void TPZDarcyAnalysis::SetLastState()
{

}

void TPZDarcyAnalysis::SetNextState()
{

}

void TPZDarcyAnalysis::Assemble()
{

}

void TPZDarcyAnalysis::AssembleLastStep()
{

}

void TPZDarcyAnalysis::AssembleResidual()
{

}

void TPZDarcyAnalysis::Residual(TPZFMatrix<STATE> &residual, int icase)
{
	TPZNonLinearAnalysis::Residual(residual, icase);
	residual = fResidualLastState + residual;
}

void TPZDarcyAnalysis::ComputeTangent(TPZFMatrix<STATE> &tangent, TPZVec<REAL> &coefs, int icase)
{
	this->SetNextState();
	TPZDarcyAnalysis::ComputeTangent(tangent, coefs, icase);
}

void TPZDarcyAnalysis::TimeForward(TPZFMatrix<STATE> &AlphasAtNplusOne, TPZFMatrix<STATE> &AlphasAtN)
{
	this->LoadSolution(AlphasAtN);
	this->AssembleLastStep();
	
	REAL tol = 1.e-6;
	int numiter = 2;
	
	this->IterativeProcess(std::cout,tol,numiter,false,false);
	
	AlphasAtNplusOne = fSolution;
}

void TPZDarcyAnalysis::InitializeFirstSolution(TPZFMatrix<STATE> &AlphasAtN, REAL &ReferencePressure)
{

}

void TPZDarcyAnalysis::Run()
{
    
    std::string dirname = PZSOURCEDIR;
    #ifdef LOG4CXX
        std::string FileName = dirname;
        FileName = dirname + "/Projects/DarcyflowAxisymmetricHdiv/";
        FileName += "DarcyFlowLog.cfg";
        InitializePZLOG(FileName);
    #endif
    
    
    //  Reading mesh
    std::string GridFileName;
    GridFileName = dirname + "/Projects/DarcyflowAxisymmetricHdiv/";
    GridFileName += "SingleLayer.dump";
//    GridFileName += "BatatacoarseQ.dump";
    
    ReadGeoMesh(GridFileName);
    RotateGeomesh(0.0*M_PI/8.0);
    this->UniformRefinement(2);
    this->PrintGeoMesh();
    
    int q = 2;
    int p = 2;
    CreateMultiphysicsMesh(q,p);
    CreateInterfaces();    // insert interfaces between bc and domain
    
    // Analysis
    bool mustOptimizeBandwidth = false;
    TPZAnalysis *an = new TPZAnalysis(fcmeshMixed,mustOptimizeBandwidth);
    TPZSkylineNSymStructMatrix skyl(fcmeshMixed);
    skyl.SetNumThreads(2);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELU);
    an->SetSolver(step);
    an->SetStructuralMatrix(skyl);
    
    this->PrintLS(an);
    
    NewtonIterations(an);
    
    this->PostProcessVTK(an);
    
}


void TPZDarcyAnalysis::CreateInterfaces()
{
    fgmesh->ResetReference();
    fcmeshMixed->LoadReferences();
    
    // Creation of interface elements
    int nel = fcmeshMixed->ElementVec().NElements();
    for(int el = 0; el < nel; el++)
    {
        TPZCompEl * compEl = fcmeshMixed->ElementVec()[el];
        if(!compEl) continue;
        int index = compEl ->Index();
        if(compEl->Dimension() == fcmeshMixed->Dimension() - 1)
        {
            TPZMultiphysicsElement * InterpEl = dynamic_cast<TPZMultiphysicsElement *>(fcmeshMixed->ElementVec()[index]);
            if(!InterpEl) continue;
            InterpEl->CreateInterfaces();
        }
    }
}


void TPZDarcyAnalysis::PrintLS(TPZAnalysis *an)
{
    an->Assemble();
    TPZAutoPointer< TPZMatrix<REAL> > matK;
    TPZFMatrix<STATE> fvec;
    matK=an->Solver().Matrix();
    fvec = an->Rhs();
    
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        matK->Print("matK = ", sout,EMathematicaInput);
        fvec.Print("fvec = ", sout,EMathematicaInput);
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
    
}

void TPZDarcyAnalysis::CreateMultiphysicsMesh(int q, int p)
{
    fmeshvec[0] = CmeshFlux(q);
    fmeshvec[1] = CmeshPressure(p);
    
    fcmeshMixed = CmeshMixed();
    
    // Transferindo para a multifisica
    TPZBuildMultiphysicsMesh::AddElements(fmeshvec, fcmeshMixed);
    TPZBuildMultiphysicsMesh::AddConnects(fmeshvec, fcmeshMixed);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(fmeshvec, fcmeshMixed);
    
    std::ofstream dumpfile("ComputationaMeshMultiphysics.txt");
    fcmeshMixed->Print(dumpfile);
    
}

void TPZDarcyAnalysis::NewtonIterations(TPZAnalysis *an)
{

    TPZFMatrix<STATE> Residual(an->Rhs().Rows(),1,0.0);
    
    an->Assemble();
    Residual = an->Rhs();
    
    TPZFMatrix<STATE> X(an->Rhs().Rows(),1,0.0);
    TPZFMatrix<STATE> DeltaX(an->Rhs().Rows(),1,0.0);
    
    STATE error=1;
    STATE normdx=1;
    int iterations=0;
    
    while (error >= fSimulationData->GetToleranceRes() && iterations <= fSimulationData->GetMaxiterations()) {
        
        an->Rhs() *= -1.0;
        an->Solve();
        DeltaX = an->Solution();
        normdx = Norm(DeltaX);
        X += DeltaX;
        
#ifdef LOG4CXX
        if(logger->isDebugEnabled())
        {
            std::stringstream sout;
            DeltaX.Print("DeltaX = ", sout,EMathematicaInput);
            X.Print("X = ", sout,EMathematicaInput);
            LOGPZ_DEBUG(logger,sout.str())
        }
#endif
        
        fcmeshMixed->LoadSolution(X);
        TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fcmeshMixed);
        an->Assemble(); // or residual
        Residual = an->Rhs();
        
        error = Norm(Residual);
        iterations++;
	
	if(error <= fSimulationData->GetToleranceRes())
	{	  
	  std::cout << "Converged with iterations:  " << iterations << std::endl;
	  std::cout << "error norm: " << error << std::endl;  
	  std::cout << "error of dx: " << normdx << std::endl;		  
	}
        
        if (iterations == fSimulationData->GetMaxiterations()) {
            std::cout << "Out max iterations " << iterations << std::endl;
            std::cout << "error norm " << error << std::endl;	    
            break;
        }
        

    }
    

}

TPZCompMesh * TPZDarcyAnalysis::CmeshMixed()
{
    int dim = 2;
    int ilayer = 0;
    int RockId = fLayers[ilayer]->GetMatIDs()[0];
    int bottomId = fLayers[ilayer]->GetMatIDs()[1];
    int rigthId = fLayers[ilayer]->GetMatIDs()[2];
    int topId = fLayers[ilayer]->GetMatIDs()[3];
    int leftId = fLayers[ilayer]->GetMatIDs()[4];
    
    const int typeFlux = 1, typePressure = 0;
    TPZFMatrix<STATE> val1(1,2,0.), val2(1,1,0.);
    
    // Malha computacional
    TPZCompMesh *cmesh = new TPZCompMesh(fgmesh);
    
    // Material medio poroso
    TPZAxiSymmetricDarcyFlow * mat = new TPZAxiSymmetricDarcyFlow(RockId);
    mat->SetReservoirData(fLayers[ilayer]);
    cmesh->InsertMaterialObject(mat);
    
//    
//    TPZAutoPointer<TPZFunction<STATE> > TimeDepFExact = new TPZDummyFunction<STATE>(PressureAnalytic);
//    mat->SetTimeDependentFunctionExact(TimeDepFExact);
    
    // Bc Bottom
    val2(0,0) = 0.0;
    TPZBndCond * bcBottom = mat->CreateBC(mat, bottomId, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(bcBottom);
    
    // Bc Right
    val2(0,0) = 10.0*1e6;
    TPZBndCond * bcRight = mat->CreateBC(mat, rigthId, typePressure, val1, val2);
    cmesh->InsertMaterialObject(bcRight);
    
    // Bc Top
    val2(0,0) = 0.0;
    TPZBndCond * bcTop = mat->CreateBC(mat, topId, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(bcTop);
    
    // Bc Left
    val2(0,0) = -0.00000001;
    TPZBndCond * bcLeft = mat->CreateBC(mat, leftId, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(bcLeft);
    

    cmesh->SetDimModel(dim);
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    cmesh->AutoBuild();
    
    
    return cmesh;
}



TPZCompMesh * TPZDarcyAnalysis::CmeshFlux(int qorder)
{

    int dim = 2;
    int ilayer = 0;
    int RockId = fLayers[ilayer]->GetMatIDs()[0];
    int bottomId = fLayers[ilayer]->GetMatIDs()[1];
    int rigthId = fLayers[ilayer]->GetMatIDs()[2];
    int topId = fLayers[ilayer]->GetMatIDs()[3];
    int leftId = fLayers[ilayer]->GetMatIDs()[4];

    const int typeFlux = 0, typePressure = 1;
    TPZFMatrix<STATE> val1(3,2,0.), val2(3,1,0.);
    
    // Malha computacional
    TPZCompMesh *cmesh = new TPZCompMesh(fgmesh);
    
    TPZAxiSymmetricDarcyFlow * mat = new TPZAxiSymmetricDarcyFlow(RockId);
    cmesh->InsertMaterialObject(mat);
    
    // Bc Bottom
    TPZBndCond * bcBottom = mat->CreateBC(mat, bottomId, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(bcBottom);
    
    // Bc Right
    TPZBndCond * bcRight = mat->CreateBC(mat, rigthId, typePressure, val1, val2);
    cmesh->InsertMaterialObject(bcRight);
    
    // Bc Top
    TPZBndCond * bcTop = mat->CreateBC(mat, topId, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(bcTop);
    
    // Bc Left
    TPZBndCond * bcLeft = mat->CreateBC(mat, leftId, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(bcLeft);
    
    // Setando Hdiv
    cmesh->SetDimModel(2);
    cmesh->SetDefaultOrder(qorder);
    cmesh->SetAllCreateFunctionsHDiv();
    
    
    cmesh->AutoBuild();

    
#ifdef DEBUG
    std::ofstream out("cmeshFlux.txt");
    cmesh->Print(out);
#endif
    
    return cmesh;
}

TPZCompMesh * TPZDarcyAnalysis::CmeshPressure(int porder)
{
    
    int dim = 2;
    int ilayer = 0;
    int RockId = fLayers[ilayer]->GetMatIDs()[0];
    int bottomId = fLayers[ilayer]->GetMatIDs()[1];
    int rigthId = fLayers[ilayer]->GetMatIDs()[2];
    int topId = fLayers[ilayer]->GetMatIDs()[3];
    int leftId = fLayers[ilayer]->GetMatIDs()[4];
    
    const int typeFlux = 0, typePressure = 1;
    TPZFMatrix<STATE> val1(3,2,0.), val2(3,1,0.);
    
    // Malha computacional
    TPZCompMesh *cmesh = new TPZCompMesh(fgmesh);
    
    TPZAxiSymmetricDarcyFlow * mat = new TPZAxiSymmetricDarcyFlow(RockId);
    cmesh->InsertMaterialObject(mat);
    
    // Bc Bottom
    TPZBndCond * bcBottom = mat->CreateBC(mat, bottomId, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(bcBottom);
    
    // Bc Right
    TPZBndCond * bcRight = mat->CreateBC(mat, rigthId, typePressure, val1, val2);
    cmesh->InsertMaterialObject(bcRight);
    
    // Bc Top
    TPZBndCond * bcTop = mat->CreateBC(mat, topId, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(bcTop);
    
    // Bc Left
    TPZBndCond * bcLeft = mat->CreateBC(mat, leftId, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(bcLeft);
    
    // Setando L2
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(porder);
    cmesh->SetAllCreateFunctionsDiscontinuous();
    cmesh->AutoBuild();
    
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
#ifdef DEBUG
    std::ofstream out("cmeshPress.txt");
    cmesh->Print(out);
#endif
    
    return cmesh;
    
}


void TPZDarcyAnalysis::ReadGeoMesh(std::string GridFileName)
{
    TPZReadGIDGrid GeometryInfo;
    GeometryInfo.SetfDimensionlessL(1.0);
    fgmesh = GeometryInfo.GeometricGIDMesh(GridFileName);
    REAL angle = 0.0*M_PI/4.0;
    RotateGeomesh(angle);
}

void TPZDarcyAnalysis::PrintGeoMesh()
{
    
#ifdef DEBUG
    //  Print Geometrical Base Mesh
    std::ofstream argument("GeometicMesh.txt");
    fgmesh->Print(argument);
    std::ofstream Dummyfile("GeometricMesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(fgmesh,Dummyfile, true);
    
#endif
}

void TPZDarcyAnalysis::RotateGeomesh(REAL CounterClockwiseAngle)
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
    
    int NumberofGeoNodes = fgmesh->NNodes();
    for (int inode = 0; inode < NumberofGeoNodes; inode++)
    {
        TPZGeoNode GeoNode = fgmesh->NodeVec()[inode];
        GeoNode.GetCoordinates(iCoords);
        // Apply rotation
        iCoordsRotated[0] = RotationMatrix(0,0)*iCoords[0]+RotationMatrix(0,1)*iCoords[1]+RotationMatrix(0,2)*iCoords[2];
        iCoordsRotated[1] = RotationMatrix(1,0)*iCoords[0]+RotationMatrix(1,1)*iCoords[1]+RotationMatrix(1,2)*iCoords[2];
        iCoordsRotated[2] = RotationMatrix(2,0)*iCoords[0]+RotationMatrix(2,1)*iCoords[1]+RotationMatrix(2,2)*iCoords[2];
        GeoNode.SetCoord(iCoordsRotated);
        fgmesh->NodeVec()[inode] = GeoNode;
    }
}

void TPZDarcyAnalysis::UniformRefinement(int nh)
{
    for ( int ref = 0; ref < nh; ref++ ){
        TPZVec<TPZGeoEl *> filhos;
        long n = fgmesh->NElements();
        for ( long i = 0; i < n; i++ ){
            TPZGeoEl * gel = fgmesh->ElementVec() [i];
            if (gel->Dimension() == 2 || gel->Dimension() == 1) gel->Divide (filhos);
        }//for i
    }//ref
}

void TPZDarcyAnalysis::PostProcessVTK(TPZAnalysis *an)
{
    const int dim = 2;
    int div = 2;
    TPZStack<std::string> scalnames, vecnames;
    std::string plotfile = "2DMixedDarcy.vtk";
    scalnames.Push("Pressure");
    vecnames.Push("Velocity");
    an->DefineGraphMesh(dim, scalnames, vecnames, plotfile);
    an->PostProcess(div);
}


