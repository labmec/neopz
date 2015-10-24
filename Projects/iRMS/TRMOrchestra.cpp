//
//  TRMOrchestra.cpp
//  PZ
//
//  Created by omar duran on 5/05/2015.
//
//

#include "TRMOrchestra.h"
#include "pzmaterial.h"
#include "TPZVecL2.h"
#include "pzl2projection.h"
#include "pzbndcond.h"
#include "TRMFlowConstants.h"
#include "pzquad.h"
#include "pzaxestools.h"


/** @brief Default constructor */
TRMOrchestra::TRMOrchestra(){
    
    /** @brief Define the global geometry being used */
    fgmesh = NULL;
    
//    /** @brief Define the mixed system analysis */
//    fFluxPressureAnalysis = NULL;
//    
//    /** @brief Define the transpor equation analysis */
//    fTransportAnalysis = NULL;
    
    /** @brief Define simulation data */
    fSimulationData = NULL;
    
}

/** @brief Default desconstructor */
TRMOrchestra::~TRMOrchestra(){
    
}

/** @brief Create a primal analysis using space odissey */
void TRMOrchestra::CreateAnalysisPrimal(){
    
    TPZManVector<int,2> dx(2,1), dy(2,1), dz(2,1);
    dx[0] = 1;
    dy[0] = 1;
    dz[0] = 1;
    
    fSpaceGenerator.CreateGeometricBoxMesh(dx, dy, dz);
//    spacegenerator.CreateGeometricReservoirMesh();
    fSpaceGenerator.PrintGeometry();
    fgmesh = fSpaceGenerator.GetGmesh();
    fSpaceGenerator.CreateH1Cmesh();
    
    TPZAutoPointer<TPZCompMesh > Cmesh = fSpaceGenerator.GetH1Cmesh();
    
    // Analysis
    bool mustOptimizeBandwidth = true;
    TPZAnalysis * AnalysisPrimal = new TPZAnalysis(Cmesh.operator->(),mustOptimizeBandwidth);
    int numofThreads = 8;
    
    TPZSkylineStructMatrix skylnsym(Cmesh.operator->());
    TPZStepSolver<STATE> step;
    skylnsym.SetNumThreads(numofThreads);
    step.SetDirect(ECholesky);
    AnalysisPrimal->SetStructuralMatrix(skylnsym);
    AnalysisPrimal->SetSolver(step);;
    AnalysisPrimal->Run();
    std::cout << "Primal dof: " << AnalysisPrimal->Rhs().Rows() << std::endl;
    
    const int dim = 3;
    int div = 1;
    TPZStack<std::string> scalnames, vecnames;
    std::string plotfile =  "../PrimalDarcy.vtk";
    scalnames.Push("Pressure");
    vecnames.Push("MinusKGradU");
    AnalysisPrimal->DefineGraphMesh(dim, scalnames, vecnames, plotfile);
    AnalysisPrimal->PostProcess(div);
    
}

/** @brief Create a dual analysis using space odissey */
void TRMOrchestra::CreateAnalysisDual(){
    int nel = 1;
    TPZManVector<int,2> dx(2,nel), dy(2,nel), dz(2,nel);
    dx[0] = 1;
    dy[0] = 1;
    dz[0] = 1;
    
//    fSpaceGenerator.CreateGeometricBoxMesh(dx, dy, dz);
    fSpaceGenerator.CreateGeometricReservoirMesh();
#ifdef PZDEBUG
    fSpaceGenerator.PrintGeometry();
#endif
    
    fSpaceGenerator.SetDefaultPOrder(2);
    
    fSpaceGenerator.CreateFluxCmesh();
    fSpaceGenerator.CreatePressureCmesh();

    
    fSpaceGenerator.CreateMixedCmesh();
    
//    ProjectExactSolution();
    
    
    fSpaceGenerator.IncreaseOrderAroundWell(3);

    fSpaceGenerator.ConfigureWellConstantPressure(0., 1000.);

    fSpaceGenerator.StaticallyCondenseEquations();
    
    // transfer the solution from the meshes to the multiphysics mesh
    TPZManVector<TPZAutoPointer<TPZCompMesh>,3 > meshvec(2);
    meshvec[0] = fSpaceGenerator.GetFluxCmesh();
    meshvec[1] = fSpaceGenerator.GetPressureMesh();
    TPZAutoPointer<TPZCompMesh > Cmesh = fSpaceGenerator.GetMixedCmesh();
    
    
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, Cmesh);
    
    
    // Analysis
    bool mustOptimizeBandwidth = true;
    fFluxPressureAnalysis.SetCompMesh(Cmesh.operator->(), mustOptimizeBandwidth);
//    TPZAnalysis * AnalysisDual = new TPZAnalysis(Cmesh.operator->(),mustOptimizeBandwidth);
    int numofThreads = 8;
#ifdef PZDEBUG
    {
        std::ofstream out("../MFCompMesh.txt");
        Cmesh->Print(out);
    }
#endif
    fFluxPressureAnalysis.Solution().Zero();
    fFluxPressureAnalysis.LoadSolution();
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, Cmesh);
    TPZFMatrix<STATE> prevsol = fFluxPressureAnalysis.Solution();
    std::cout << "Total dof: " << prevsol.Rows() << std::endl;
    
    std::map<REAL,STATE> RatebyPosition;
    STATE TotalFlux = 0.;
     ComputeProductionRate(RatebyPosition, TotalFlux);

    TPZSkylineStructMatrix strmat(Cmesh.operator->());
//    TPZSkylineNSymStructMatrix strmat(Cmesh.operator->());
    TPZStepSolver<STATE> step;
    strmat.SetNumThreads(numofThreads);
    step.SetDirect(ELDLt);
    fFluxPressureAnalysis.SetStructuralMatrix(strmat);
    fFluxPressureAnalysis.SetSolver(step);
    fFluxPressureAnalysis.Run();
    std::cout << "Rhs norm " << Norm(fFluxPressureAnalysis.Rhs()) << std::endl;
    prevsol -= fFluxPressureAnalysis.Solution();
    Cmesh->LoadSolution(prevsol);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, Cmesh);

    RatebyPosition.clear();
    TotalFlux = 0.;
    ComputeProductionRate(RatebyPosition, TotalFlux);    
    std::cout << "Total flux " << TotalFlux << std::endl;
    for (std::map<REAL,STATE>::iterator it = RatebyPosition.begin(); it != RatebyPosition.end(); it++) {
        std::cout << "Y_position " << it->first << " rate " << it->second << std::endl;
    }
    
    if (0)
    {
        fFluxPressureAnalysis.Run();
        std::cout << "Rhs norm " << Norm(fFluxPressureAnalysis.Rhs()) << std::endl;
        prevsol -= fFluxPressureAnalysis.Solution();
        Cmesh->LoadSolution(prevsol);
        TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, Cmesh);
        fFluxPressureAnalysis.AssembleResidual();
        std::cout << "Rhs norm " << Norm(fFluxPressureAnalysis.Rhs()) << std::endl;
    }
    const int dim = 3;
    int div = 0;
    TPZStack<std::string> scalnames, vecnames;
    std::string plotfile =  "../DualDarcy.vtk";
    scalnames.Push("WeightedPressure");
    scalnames.Push("DivOfBulkVeclocity");
    scalnames.Push("POrder");
    vecnames.Push("BulkVelocity");
    fFluxPressureAnalysis.DefineGraphMesh(dim, scalnames, vecnames, plotfile);
    fFluxPressureAnalysis.PostProcess(div);
    
}

/** @brief Create computational meshes using space odissey */
void TRMOrchestra::CreateCompMeshes(TRMRawData &rawdata){
    
}

/// Transfer the flux solution to the saturation mesh
void TRMOrchestra::TransferToSaturationMesh()
{
    //
}

/** @brief Project an exact solution */
void TRMOrchestra::ProjectExactSolution()
{
    TPZAutoPointer<TPZCompMesh> mesh = fSpaceGenerator.GetFluxCmesh();
    if (!mesh) {
        DebugStop();
    }
    // copiar os materiais
    std::map<int, TPZMaterial *> matmap;// = mesh->MaterialVec();
    std::map<int, TPZMaterial *>::iterator it;
    for (it= mesh->MaterialVec().begin(); it != mesh->MaterialVec().end(); it++) {
        it->second->Clone(matmap);
    }
    
    // put L2 projection material
    TPZVecL2 *vecmat = new TPZVecL2(_ReservMatId);
    TPZAutoPointer<TPZFunction<STATE> > force = new TPZDummyFunction<STATE>(TRMOrchestra::ExactFlux);
    vecmat->SetForcingFunction(force);
    for (it = mesh->MaterialVec().begin(); it != mesh->MaterialVec().end(); it++) {
        delete it->second;
    }
    mesh->MaterialVec().clear();
    mesh->InsertMaterialObject(vecmat);
    {
        TPZAnalysis an(mesh,false);
        TPZSkylineStructMatrix strmat(mesh);
        std::set<int> matids;
        matids.insert(_ReservMatId);
        strmat.SetMaterialIds(matids);
        an.SetStructuralMatrix(strmat);
        TPZStepSolver<STATE> step;
        step.SetDirect(ELDLt);
        an.SetSolver(step);
        an.Run();
        TPZStack<std::string> scalnames,vecnames;
        vecnames.Push("state");
        an.DefineGraphMesh(3, scalnames, vecnames, "../ProjectFlux.vtk");
        an.PostProcess(0);
    }
    delete vecmat;
    mesh->MaterialVec() = matmap;
    
    mesh = fSpaceGenerator.GetPressureMesh();
    matmap.clear();// = mesh->MaterialVec();
    for (it= mesh->MaterialVec().begin(); it != mesh->MaterialVec().end(); it++) {
        it->second->Clone(matmap);
    }
    int nstate = 1;
    TPZManVector<STATE,1> sol(1,1.);
    int dim = 3;
    TPZL2Projection *l2proj = new TPZL2Projection(_ReservMatId,dim,nstate,sol);
    TPZAutoPointer<TPZFunction<STATE> > force2 = new TPZDummyFunction<STATE>(TRMOrchestra::ExactPressure);
    l2proj->SetForcingFunction(force2);

    for (it = mesh->MaterialVec().begin(); it != mesh->MaterialVec().end(); it++) {
        delete it->second;
    }
    mesh->MaterialVec().clear();
    mesh->InsertMaterialObject(l2proj);
    {
        TPZAnalysis an(mesh,false);
        TPZSkylineStructMatrix strmat(mesh);
        std::set<int> matids;
        matids.insert(_ReservMatId);
        strmat.SetMaterialIds(matids);
        an.SetStructuralMatrix(strmat);
        TPZStepSolver<STATE> step;
        step.SetDirect(ELDLt);
        an.SetSolver(step);
        an.Run();
//        mesh->Print(std::cout);
        TPZStack<std::string> scalnames,vecnames;
        scalnames.Push("state");
        an.DefineGraphMesh(3, scalnames, vecnames, "../ProjectPressure.vtk");
        an.PostProcess(0);
    }
    delete l2proj;
    mesh->MaterialVec() = matmap;
    
    mesh = fSpaceGenerator.GetMixedCmesh();
    for (it = mesh->MaterialVec().begin(); it != mesh->MaterialVec().end(); it++) {
        TPZBndCond * bnd = dynamic_cast<TPZBndCond *>(it->second);
        if (bnd)
        {
            bnd->SetForcingFunction(0,force2);
        }
    }

}

/** @brief exact pressure */
void TRMOrchestra::ExactPressure(const TPZVec<REAL> &x, TPZVec<STATE> &pressure)
{
    pressure[0] = 0.;//x[2];
}

/** @brief exact pressure */
void TRMOrchestra::ExactFlux(const TPZVec<REAL> &x, TPZVec<STATE> &flux)
{
    flux[0] = 0.;
    flux[1] = 0.;
    flux[2] = 0.;//-1.;
    
}

/** @brief Compute the production rate of the reservoir */
void TRMOrchestra::ComputeProductionRate(std::map<REAL,STATE> &RatebyPosition, STATE &TotalIntegral)
{
    fSpaceGenerator.GetGmesh()->ResetReference();
    fSpaceGenerator.GetFluxCmesh()->LoadReferences();
    TPZAutoPointer<TPZGeoMesh> gmesh = fSpaceGenerator.GetGmesh();
    TPZInt1d intrule(14);
    int np = intrule.NPoints();
    TotalIntegral = 0.;
    long nel = gmesh->NElements();
    for (long el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (gel->MaterialId() != _WellMatId3D) {
            continue;
        }
//        std::cout << "Accumulating for well3d = " << gel->Index() << std::endl;
        TPZStack<long> faceElements;
        for (int side=21; side < 25 ; side++) {
            TPZGeoElSide gelside(gel,side);
            TPZGeoElSide neighbour = gelside.Neighbour();
            while (neighbour.Element()->MaterialId() != _Well3DReservoirFaces && neighbour != gelside) {
                neighbour = neighbour.Neighbour();
            }
            if (neighbour == gelside) {
                DebugStop();
            }
            faceElements.Push(neighbour.Element()->Index());
        }
        for (int i=0; i<np; i++) {
            REAL yIntegral = 0.;
            REAL lineIntegral = 0.;
            STATE flux1DIntegral = 0.;
            TPZManVector<REAL,2> posi(1);
            TPZStack<REAL> yvals;
            REAL wi;
            intrule.Point(i,posi,wi);
            for (int iface=0; iface< faceElements.size(); iface++) {
                for (int j=0; j<np; j++) {
                    TPZGeoEl *gelface = gmesh->Element(faceElements[iface]);
                    TPZCompEl *celface = gelface->Reference();
                    TPZFNMatrix<9,REAL> jac(2,2),jacinv(2,2),axes(2,3),gradx(3,2);
                    REAL detjac,wj;
                    TPZManVector<REAL,3> posj(1),intpoint(2),xco(3);
                    intrule.Point(j, posj, wj);
                    intpoint[1] = posi[0];
                    intpoint[0] = posj[0];
                    gelface->Jacobian(intpoint, jac, axes, detjac, jacinv);
                    gelface->X(intpoint, xco);
//                    std::cout << "el index " << gelface->Index() << " intpoint " << intpoint << " xco " << xco << std::endl;
                    TPZAxesTools<REAL>::ComputeGradX(jac, axes, gradx);
                    int xdir = 0;
                    REAL detjac1d = sqrt(gradx(0,xdir)*gradx(0,xdir)+gradx(1,xdir)*gradx(1,xdir)+gradx(2,xdir)*gradx(2,xdir));
                    lineIntegral += detjac1d*wj;
                    yIntegral += detjac1d*wj*xco[1];
                    yvals.Push(xco[1]);
                    TPZManVector<STATE,1> sol(1);
                    celface->Solution(intpoint, 0, sol);
                    flux1DIntegral += detjac1d*wj*sol[0];
                    TotalIntegral += fabs(detjac)*wi*wj*sol[0];
                }
            }
//            std::cout << "Y values " << yvals << std::endl;
            REAL averageY = yIntegral/lineIntegral;
            RatebyPosition[averageY] = flux1DIntegral;
        }
    }
}
