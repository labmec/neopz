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
#include "tpzgeoelrefpattern.h"
#include "TPZGeoLinear.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"

#include "TPZVTKGeoMesh.h"
#include "TPZAxiSymmetricDarcyFlow.h"
#include "TPZAxiSymmetricDarcyFlowH1.h"
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
    fcmeshdarcy=NULL;
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
//    GridFileName += "MixLayer.dump";
//    GridFileName += "BatatacoarseQ.dump";
    
    if(fLayers[0]->GetIsGIDGeometry())
    {
        ReadGeoMesh(GridFileName);
    }
    else
    {
        CreatedGeoMesh();
    }
    

    RotateGeomesh(0.0*M_PI/8.0);
    this->UniformRefinement(1);
    this->PrintGeoMesh();
    
    
    
    int q = 2;
    int p = 2;
    
    if (fSimulationData->GetIsH1approx())
    {
        CmeshH1(p);
    }
    else
    {
        CreateMultiphysicsMesh(q,p);
        CreateInterfaces();    // insert interfaces between bc and domain
    }

    
    // Analysis
    bool mustOptimizeBandwidth = true;
    TPZAnalysis *an = new TPZAnalysis(fcmeshdarcy,mustOptimizeBandwidth);
    TPZSkylineNSymStructMatrix skyl(fcmeshdarcy);
    skyl.SetNumThreads(4);
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
    fcmeshdarcy->LoadReferences();
    
    // Creation of interface elements
    int nel = fcmeshdarcy->ElementVec().NElements();
    for(int el = 0; el < nel; el++)
    {
        TPZCompEl * compEl = fcmeshdarcy->ElementVec()[el];
        if(!compEl) continue;
        int index = compEl ->Index();
        if(compEl->Dimension() == fcmeshdarcy->Dimension() - 1)
        {
            TPZMultiphysicsElement * InterpEl = dynamic_cast<TPZMultiphysicsElement *>(fcmeshdarcy->ElementVec()[index]);
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
    
    fcmeshdarcy = CmeshMixed();
    
    // Transferindo para a multifisica
    TPZBuildMultiphysicsMesh::AddElements(fmeshvec, fcmeshdarcy);
    TPZBuildMultiphysicsMesh::AddConnects(fmeshvec, fcmeshdarcy);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(fmeshvec, fcmeshdarcy);
    
    std::ofstream dumpfile("ComputationaMeshMultiphysics.txt");
    fcmeshdarcy->Print(dumpfile);
    
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
        
        fcmeshdarcy->LoadSolution(X);
        if (!fSimulationData->GetIsH1approx())
        {
            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fcmeshdarcy);
        }

        
        an->Assemble(); // or residual
        Residual = an->Rhs();
        
#ifdef LOG4CXX
        if(logger->isDebugEnabled())
        {
            std::stringstream sout;
            Residual.Print("Residual = ", sout,EMathematicaInput);
            LOGPZ_DEBUG(logger,sout.str())
        }
#endif
        
        error = Norm(Residual);
        iterations++;
	
	if(error <= fSimulationData->GetToleranceRes() || normdx <= fSimulationData->GetToleranceDX())
	{	  
	  std::cout << "Converged with iterations:  " << iterations << std::endl;
	  std::cout << "error norm: " << error << std::endl;  
	  std::cout << "error of dx: " << normdx << std::endl;
      break;
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
    val2(0,0) = -0.00002;
    TPZBndCond * bcLeft = mat->CreateBC(mat, leftId, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(bcLeft);
    

    cmesh->SetDimModel(dim);
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    cmesh->AutoBuild();
    
    
    return cmesh;
}

void TPZDarcyAnalysis::CmeshH1(int porder)
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
    TPZAxiSymmetricDarcyFlowH1 * mat = new TPZAxiSymmetricDarcyFlowH1(RockId);
    mat->SetReservoirData(fLayers[ilayer]);
    cmesh->InsertMaterialObject(mat);
    
    
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
    val2(0,0) = -0.00002;
    TPZBndCond * bcLeft = mat->CreateBC(mat, leftId, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(bcLeft);
    
    
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(porder);
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->AutoBuild();
    
#ifdef DEBUG
    std::ofstream out("cmeshPressureH1.txt");
    cmesh->Print(out);
#endif
    
    fcmeshdarcy = cmesh;

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
    GeometryInfo.SetfDimensionlessL(0.001);
    fgmesh = GeometryInfo.GeometricGIDMesh(GridFileName);
    REAL angle = 0.0*M_PI/4.0;
    RotateGeomesh(angle);
}

void TPZDarcyAnalysis::CreatedGeoMesh()
{
    long Qnodes = 4;
    int ilayer = 0;

    TPZGeoMesh *gmesh= new TPZGeoMesh;
    
    gmesh->SetMaxNodeId(Qnodes-1);
    gmesh->NodeVec().Resize(Qnodes);
    TPZVec<TPZGeoNode> Node(Qnodes);
    
    TPZVec <long> TopolQuad(4);
    TPZVec <long> TopolLine(2);
    REAL r     = fLayers[ilayer]->Layerr();
    REAL rw    = fLayers[ilayer]->Layerrw();
    REAL h     = fLayers[ilayer]->Layerh();
    REAL top   = fLayers[ilayer]->LayerTop();
    
    int RockId = fLayers[ilayer]->GetMatIDs()[0];
    int bottomId = fLayers[ilayer]->GetMatIDs()[1];
    int rigthId = fLayers[ilayer]->GetMatIDs()[2];
    int topId = fLayers[ilayer]->GetMatIDs()[3];
    int leftId = fLayers[ilayer]->GetMatIDs()[4];
    
    // Nodes
    long id = 0;
    
    Node[id].SetNodeId(id);
    Node[id].SetCoord(0 ,  rw);         //coord r
    Node[id].SetCoord(1 , top - h);     //coord z
    gmesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    Node[id].SetCoord(0 ,  rw + r);         //coord r
    Node[id].SetCoord(1 , top - h);     //coord z
    gmesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    Node[id].SetCoord(0 ,  rw + r);         //coord r
    Node[id].SetCoord(1 ,  top);     //coord z
    gmesh->NodeVec()[id] = Node[id];
    id++;
    
    Node[id].SetNodeId(id);
    Node[id].SetCoord(0 , rw);         //coord r
    Node[id].SetCoord(1 , top);     //coord z
    gmesh->NodeVec()[id] = Node[id];
    id++;
    
    
//  Geometric Elements
    int elid = 0;
    
    TopolLine[0] = 0;
    TopolLine[1] = 1;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elid,TopolLine,bottomId,*gmesh);
    id++;
    
    TopolLine[0] = 1;
    TopolLine[1] = 2;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elid,TopolLine,rigthId,*gmesh);
    id++;
    
    TopolLine[0] = 2;
    TopolLine[1] = 3;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elid,TopolLine,topId,*gmesh);
    id++;
    
    TopolLine[0] = 0;
    TopolLine[1] = 3;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elid,TopolLine,leftId,*gmesh);
    id++;
    
    
    TopolQuad[0] = 0;
    TopolQuad[1] = 1;
    TopolQuad[2] = 2;
    TopolQuad[3] = 3;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (elid,TopolQuad,RockId,*gmesh);
    
    
    gmesh->BuildConnectivity();
    fgmesh = gmesh;

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
    std::string plotfile;
    if (fSimulationData->GetIsH1approx()) {
        plotfile = "2DH1Darcy.vtk";
    }
    else{
        plotfile = "2DMixedDarcy.vtk";
    }

    scalnames.Push("Pressure");
    scalnames.Push("Density");
    scalnames.Push("Porosity");
    scalnames.Push("DivofVeclocity");    
    vecnames.Push("Velocity");
    an->DefineGraphMesh(dim, scalnames, vecnames, plotfile);
    an->PostProcess(div);
}


