//
//  TPZDarcyAnalysis.cpp
//  PZ
//
//  Created by Nathan Shauer and Omar Duran on 9/8/14.
//
//

#include "pzlog.h"
#include "TPZDarcyAnalysis.h"
#include "TPZMatDarcy2dhdiv.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "TPZVTKGeoMesh.h"
#include "TPZReadGIDGrid.h"
#include "pzbndcond.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzanalysis.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzstepsolver.h"
#include "TPZCompElDisc.h"


#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.multiphase"));
#endif

TPZDarcyAnalysis::TPZDarcyAnalysis(TPZAutoPointer<TPZFracData> Data)
{
    fData = Data;
    fmeshvec.Resize(2);
}


TPZDarcyAnalysis::~TPZDarcyAnalysis()
{
    
}

void TPZDarcyAnalysis::Run()
{
    // Parametros
    const int nel = 0;
    
    // Malha geometrica
    fgmesh = CreateGMesh(nel);
    
    fmeshvec[0] = CreateCMeshFluxHdiv();
    fmeshvec[1] = CreateCMeshPressureL2();
    fcmeshMixed = CreateCMeshMixed();
    
    // Transferindo para a multifisica
	TPZBuildMultiphysicsMesh::AddElements(fmeshvec, fcmeshMixed);
	TPZBuildMultiphysicsMesh::AddConnects(fmeshvec, fcmeshMixed);
	TPZBuildMultiphysicsMesh::TransferFromMeshes(fmeshvec, fcmeshMixed);
    
    // Create Interfaces
    CreateInterfaces(fcmeshMixed);
    
    // Analysis
    bool mustOptimizeBandwidth = true;
    TPZAnalysis *an = new TPZAnalysis(fcmeshMixed,mustOptimizeBandwidth);
    TPZSkylineNSymStructMatrix skyl(fcmeshMixed);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELU);
    an->SetSolver(step);
    an->SetStructuralMatrix(skyl);
    
    SolveSistTransient(an);
    
    delete an;
    

}

TPZGeoMesh * TPZDarcyAnalysis::CreateGMesh(const int nel)
{
    std::string dirname = PZSOURCEDIR;
    std::string FileName = dirname;
    
    //  Reading mesh
    std::string GridFileName;
    GridFileName = dirname + "/Projects/OilWaterSystem/";
    GridFileName += "OilWaterSystemUnit.dump";
    REAL angle = 0.0*M_PI/4.0;
    
    TPZReadGIDGrid GeometryInfo;
    GeometryInfo.SetfDimensionlessL(100.0);
    TPZGeoMesh * gmesh = GeometryInfo.GeometricGIDMesh(GridFileName);
    RotateGeomesh(gmesh, angle);
    
    UniformRefinement(gmesh, nel);

    
#ifdef DEBUG
    //  Print Geometrical Base Mesh
    std::ofstream argument("GeometicMesh.txt");
    gmesh->Print(argument);
    std::ofstream Dummyfile("GeometricMesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh,Dummyfile, true);
#endif
    
    return gmesh;
}

TPZCompMesh * TPZDarcyAnalysis::CreateCMeshFluxHdiv()
{
    const int matId = 1, bcBottomId = 2, bcRightId = 3, bcTopId = 4, bcLeftId = 5;
    const int typeFlux = 0, typePressure = 1;
    const int fluxorder = fData->PorderFlow();
    TPZFMatrix<STATE> val1(3,2,0.), val2(3,1,0.);
    
    // Malha computacional
    TPZCompMesh *cmesh = new TPZCompMesh(fgmesh);
    
    // Material da fratura
    TPZMatDarcy2dhdiv *mat = new TPZMatDarcy2dhdiv(matId);
    cmesh->InsertMaterialObject(mat);
    
    // Bc Bottom
    TPZBndCond * bcBottom = mat->CreateBC(mat, bcBottomId, typeFlux, val1, val2);
 	cmesh->InsertMaterialObject(bcBottom);
    
    // Bc Right
    TPZBndCond * bcRight = mat->CreateBC(mat, bcRightId, typePressure, val1, val2);
 	cmesh->InsertMaterialObject(bcRight);
    
    // Bc Top
    TPZBndCond * bcTop = mat->CreateBC(mat, bcTopId, typeFlux, val1, val2);
 	cmesh->InsertMaterialObject(bcTop);
    
    // Bc Left
    TPZBndCond * bcLeft = mat->CreateBC(mat, bcLeftId, typeFlux, val1, val2);
 	cmesh->InsertMaterialObject(bcLeft);
    
    // Setando Hdiv
    cmesh->SetDimModel(2);
    cmesh->SetDefaultOrder(fluxorder);
    cmesh->SetAllCreateFunctionsHDiv();
    cmesh->AutoBuild();
    
    return cmesh;
}

TPZCompMesh * TPZDarcyAnalysis::CreateCMeshPressureL2()
{

    const int matId = 1, bcBottomId = 2, bcRightId = 3, bcTopId = 4, bcLeftId = 5;
    const int typeFlux = 0, typePressure = 1;
    const int pressureorder = fData->PorderPressure();
    TPZFMatrix<STATE> val1(3,2,0.), val2(3,1,0.);
    
    // Malha computacional
    TPZCompMesh *cmesh = new TPZCompMesh(fgmesh);
    
    // Material da fratura
    TPZMatDarcy2dhdiv *mat = new TPZMatDarcy2dhdiv(matId);
    cmesh->InsertMaterialObject(mat);
    
    // Bc Bottom
    TPZBndCond * bcBottom = mat->CreateBC(mat, bcBottomId, typeFlux, val1, val2);
 	cmesh->InsertMaterialObject(bcBottom);
    
    // Bc Right
    TPZBndCond * bcRight = mat->CreateBC(mat, bcRightId, typePressure, val1, val2);
 	cmesh->InsertMaterialObject(bcRight);
    
    // Bc Top
    TPZBndCond * bcTop = mat->CreateBC(mat, bcTopId, typeFlux, val1, val2);
 	cmesh->InsertMaterialObject(bcTop);
    
    // Bc Left
    TPZBndCond * bcLeft = mat->CreateBC(mat, bcLeftId, typeFlux, val1, val2);
 	cmesh->InsertMaterialObject(bcLeft);
    
    // Setando L2
    cmesh->SetDimModel(2);
    cmesh->SetDefaultOrder(pressureorder);
    cmesh->SetAllCreateFunctionsDiscontinuous();
    cmesh->AutoBuild();
//    cmesh->AdjustBoundaryElements();
//    cmesh->CleanUpUnconnectedNodes();
    
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
    
    return cmesh;
    
}

TPZCompMesh * TPZDarcyAnalysis::CreateCMeshMixed()
{
    // Definicao de ids e tipos
    const int matId = 1, bcBottomId = 2, bcRightId = 3, bcTopId = 4, bcLeftId = 5;
    const int typeFlux = 0, typePressure = 1;
    TPZFMatrix<STATE> val1(3,2,0.), val2(3,1,0.);
    
    // Malha computacional
    TPZCompMesh *cmesh = new TPZCompMesh(fgmesh);
    
    // Material medio poroso
    TPZMatDarcy2dhdiv *mat = new TPZMatDarcy2dhdiv(matId);
    mat->SetSimulationData(fData);
    cmesh->InsertMaterialObject(mat);
    
    // Bc Bottom
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    val2(2,0) = 0.0;
    TPZBndCond * bcBottom = mat->CreateBC(mat, bcBottomId, typeFlux, val1, val2);
 	cmesh->InsertMaterialObject(bcBottom);
    
    // Bc Right
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    val2(2,0) = 0.0;
    TPZBndCond * bcRight = mat->CreateBC(mat, bcRightId, typePressure, val1, val2);
 	cmesh->InsertMaterialObject(bcRight);
    
    // Bc Top
    val2(0,0) = 0.0;
    val2(1,0) = 0.0;
    val2(2,0) = 0.0;
    TPZBndCond * bcTop = mat->CreateBC(mat, bcTopId, typeFlux, val1, val2);
 	cmesh->InsertMaterialObject(bcTop);
    
    // Bc Left
    val2(0,0) = 1.0;
    val2(1,0) = 0.0;
    val2(2,0) = 0.0;
    TPZBndCond * bcLeft = mat->CreateBC(mat, bcLeftId, typeFlux, val1, val2);
 	cmesh->InsertMaterialObject(bcLeft);
    
    // Setando Multifisico
    cmesh->SetDimModel(2);
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    cmesh->AutoBuild();
    
    return cmesh;
}

void TPZDarcyAnalysis::CreateInterfaces(TPZCompMesh *cmesh)
{
    cmesh->Reference()->ResetReference();
    cmesh->LoadReferences();
    
    // Creation of interface elements
    int nel = cmesh->ElementVec().NElements();
    for(int el = 0; el < nel; el++)
    {
        TPZCompEl * compEl = cmesh->ElementVec()[el];
        if(!compEl) continue;
        int index = compEl ->Index();
        if(compEl->Dimension() == cmesh->Dimension())
        {
            TPZMultiphysicsElement * InterpEl = dynamic_cast<TPZMultiphysicsElement *>(cmesh->ElementVec()[index]);
            if(!InterpEl) continue;
            InterpEl->CreateInterfaces();
        }
    }
}

void TPZDarcyAnalysis::IterativeProcess(TPZAnalysis *an, std::ostream &out)
{
	int iter = 0;
	REAL error = 1.e10;
    const REAL tol = 1.e-4;
    const int numiter = 50;
    
    fData->SetCurrentState();
	int numeq = an->Mesh()->NEquations();
	
	TPZFMatrix<STATE> prevsol(an->Solution());
    TPZFMatrix<STATE> SoliterK(prevsol);
	if(prevsol.Rows() != numeq) prevsol.Redim(numeq,1);
    
    an->Assemble();
    an->Rhs() += fLastStepRhs;
    an->Rhs() *= -1.0;
    
    TPZAutoPointer< TPZMatrix<REAL> > matK;
    
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        matK=an->Solver().Matrix();
        matK->Print("matK = ", sout,EMathematicaInput);
        an->Rhs().Print("Rhs = ", sout,EMathematicaInput);
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
    
	while(error > tol && iter < numiter) {
		
		an->Solve(); // o an->Solution() eh o deltaU aqui
        SoliterK = prevsol + an->Solution();
		REAL normDeltaSol = Norm(an->Solution());
        
#ifdef LOG4CXX
        if(logger->isDebugEnabled())
        {
            std::stringstream sout;
            matK=an->Solver().Matrix();
            matK->Print("matK = ", sout,EMathematicaInput);
            an->Solution().Print("DeltaX = ", sout,EMathematicaInput);
            SoliterK.Print("Xk = ", sout,EMathematicaInput);
            LOGPZ_DEBUG(logger,sout.str())
        }
#endif
        
        an->LoadSolution(SoliterK); // Aqui o an->Solution() eh o U
        TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fcmeshMixed);
        an->Assemble();

#ifdef LOG4CXX
        if(logger->isDebugEnabled())
        {
            std::stringstream sout;
            fLastStepRhs.Print("ResAtn = ", sout,EMathematicaInput);
            an->Rhs().Print("Respone = ", sout,EMathematicaInput);            
            LOGPZ_DEBUG(logger,sout.str())
        }
#endif
        
        an->Rhs() += fLastStepRhs;
        
#ifdef LOG4CXX
        if(logger->isDebugEnabled())
        {
            std::stringstream sout;
            an->Rhs().Print("Res = ", sout,EMathematicaInput);
            LOGPZ_DEBUG(logger,sout.str())
        }
#endif
        
        double NormResLambda = Norm(an->Rhs());
		double norm = normDeltaSol;
		out << "Iteracao n : " << (iter+1) << " : normas |Delta(Un)| e |Delta(rhs)| : " << normDeltaSol << " / " << NormResLambda << std::endl;
        
		if(norm < tol) {
			out << "\nTolerancia atingida na iteracao : " << (iter+1) << std::endl;
			out << "\n\nNorma da solucao |Delta(Un)|  : " << norm << std::endl << std::endl;
            
		}
        else if( (norm - error) > 1.e-9 ) {
            out << "\nDivergent Method\n";
        }
        
		error = norm;
		iter++;
        prevsol = SoliterK;
		out.flush();
	}
    
    if (error > tol) {
        DebugStop(); // Metodo nao convergiu!!
    }
    
}

void TPZDarcyAnalysis::AssembleLastStep(TPZAnalysis *an)
{
    fData->SetLastState();
    an->Assemble();
    fLastStepRhs = an->Rhs();
}

void TPZDarcyAnalysis::SolveSistTransient(TPZAnalysis *an)
{
    
    bool mustStop = false;
    while (mustStop == false) {
        
        AssembleLastStep(an);
        
#ifdef LOG4CXX
        if(logger->isDebugEnabled())
        {
            std::stringstream sout;
            an->Rhs().Print("ResAtn = ", sout,EMathematicaInput);
            LOGPZ_DEBUG(logger,sout.str())
        }
#endif
        
        IterativeProcess(an, std::cout);
        
        fData->SetNextTime();
        
        const int dim = 2;
        int div =2;
        TPZStack<std::string> scalnames, vecnames;
        std::string plotfile = "2DMixedDarcy.vtk";
        scalnames.Push("Pressure");
        vecnames.Push("BulkVelocity");
        an->DefineGraphMesh(dim, scalnames, vecnames, plotfile);
        an->PostProcess(div,dim);
        
        REAL peteleco = 1.E-8;
        if( fData->Time() > (fData->TotalTime() - peteleco) )
        {
            mustStop = true;
        }
    }
}

void TPZDarcyAnalysis::UniformRefinement(TPZGeoMesh *gMesh, int nh)
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

void TPZDarcyAnalysis::RotateGeomesh(TPZGeoMesh *gmesh, REAL CounterClockwiseAngle)
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