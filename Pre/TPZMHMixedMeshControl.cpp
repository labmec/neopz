//
//  TPZMHMixedMeshControl.cpp
//  PZ
//
//  Created by Philippe Devloo on 09/10/16.
//
//

#include "TPZMHMixedMeshControl.h"
#include "TPZVecL2.h"
#include "pzbndcond.h"
#include "TPZNullMaterial.h"
#include "TPZMatLaplacian.h"
#include "TPZLagrangeMultiplier.h"


#include <iostream>
#include <sstream>
#include <iterator>
#include <numeric>
#include <algorithm>

#include "pzsubcmesh.h"

#include "pzbuildmultiphysicsmesh.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZCompMeshTools.h"
#include "pzinterpolationspace.h"
#include "pzintel.h"
#include "TPZInterfaceEl.h"
#include "TPZMultiphysicsInterfaceEl.h"

#include "TPZVTKGeoMesh.h"
#include "TPZNullMaterial.h"

#include "TPZGeoElSideAncestors.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mhmixedmeshcontrol"));
#endif

/*
TPZMHMixedMeshControl::TPZMHMixedMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh, std::set<int64_t> &coarseindices) : TPZMHMeshControl(gmesh,coarseindices)
{
    fFluxMesh = new TPZCompMesh(gmesh);
    fFluxMesh->SetDimModel(gmesh->Dimension());
    fCMesh = new TPZCompMesh(gmesh);
    fCMesh->SetDimModel(gmesh->Dimension());
}
*/

TPZMHMixedMeshControl::TPZMHMixedMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh, TPZVec<int64_t> &coarseindices) : TPZMHMeshControl(gmesh,coarseindices)
{
    fFluxMesh->SetName("FluxMesh");

    fRotationMesh = new TPZCompMesh(gmesh);
    fRotationMesh->SetDimModel(gmesh->Dimension());
    fDistributedFlux = new TPZCompMesh(gmesh);
    fDistributedFlux->SetDimModel(gmesh->Dimension());
    fAverageSolution = new TPZCompMesh(gmesh);
    fAverageSolution->SetDimModel(gmesh->Dimension());
    fCMesh = new TPZMultiphysicsCompMesh(gmesh);
    fCMesh->SetDimModel(gmesh->Dimension());
}


TPZMHMixedMeshControl::TPZMHMixedMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh) : TPZMHMeshControl(gmesh)
{
    fFluxMesh->SetName("FluxMesh");
    fRotationMesh = new TPZCompMesh(gmesh);
    fRotationMesh->SetDimModel(gmesh->Dimension());
    fDistributedFlux = new TPZCompMesh(gmesh);
    fDistributedFlux->SetDimModel(gmesh->Dimension());
    fAverageSolution = new TPZCompMesh(gmesh);
    fAverageSolution->SetDimModel(gmesh->Dimension());
    fCMesh = new TPZMultiphysicsCompMesh(gmesh);
    fCMesh->SetDimModel(gmesh->Dimension());
}

TPZMHMixedMeshControl::~TPZMHMixedMeshControl()
{
    if(fCMesh)
    {
        fCMesh->CleanUp();
    }
    if(fFluxMesh)
    {
        int64_t nc = fFluxMesh->NConnects();
        for (int64_t ic = 0; ic<nc; ic++) {
            TPZConnect &c = fFluxMesh->ConnectVec()[ic];
            if (c.HasDependency()) {
                c.RemoveDepend();
            }
        }
        fFluxMesh->CleanUp();
    }
    if (fPressureFineMesh)
    {
        DeletePressureElements();
        {
            int64_t nc = fPressureFineMesh->NConnects();
            for (int64_t ic = 0; ic<nc; ic++) {
                TPZConnect &c = fPressureFineMesh->ConnectVec()[ic];
                if (c.HasDependency()) {
                    c.RemoveDepend();
                }
            }
        }
        fPressureFineMesh->CleanUp();
    }
}

/// Insert Boundary condition objects that do not perform any actual computation
void TPZMHMixedMeshControl::InsertPeriferalMaterialObjects()
{
#ifdef PZDEBUG
    int matid = *fMaterialIds.begin();
    TPZMaterial *mat = fCMesh->FindMaterial(matid);
    if (!mat) {
        DebugStop();
    }
    int nstate = mat->NStateVariables();
    if(nstate != fNState) DebugStop();
#endif
    if(fCMesh->MaterialVec().find(fSkeletonMatId) != fCMesh->MaterialVec().end())
    {
        std::cout << "Peripheral material inserted twice " << __PRETTY_FUNCTION__ << std::endl;
        DebugStop();
        return;
    }
    int dim = fGMesh->Dimension();
    {
        TPZNullMaterial *nullmat = new TPZNullMaterial(fSkeletonMatId);
        nullmat->SetDimension(dim-1);
        nullmat->SetNStateVariables(fNState);
        fCMesh->InsertMaterialObject(nullmat);
    }
    {
        TPZNullMaterial *nullmat = new TPZNullMaterial(fWrapMatId);
        nullmat->SetDimension(dim-1);
        nullmat->SetNStateVariables(fNState);
        fCMesh->InsertMaterialObject(nullmat);
    }
    if(fHybridizeSkeleton)
    {
        {
            TPZNullMaterial *nullmat = new TPZNullMaterial(fSecondSkeletonMatId);
            nullmat->SetDimension(dim-1);
            nullmat->SetNStateVariables(fNState);
            fCMesh->InsertMaterialObject(nullmat);
        }
        {
            TPZNullMaterial *nullmat = new TPZNullMaterial(fPressureSkeletonMatId);
            nullmat->SetDimension(dim-1);
            nullmat->SetNStateVariables(fNState);
            fCMesh->InsertMaterialObject(nullmat);
        }
    }
    InsertPeriferalHdivMaterialObjects();
    InsertPeriferalPressureMaterialObjects();
    if(fNState > 1){
        InsertPeriferalRotationMaterialObjects();
    }
    {
        // insert the material in the constant pressure mesh and lagrange mesh
        int nStateVariables = 1;
        if (fProblemType == EScalar) {
            nStateVariables = 1;
        }
        else if(fProblemType == EElasticity2D)
        {
            nStateVariables = 3;
        }
        else if(fProblemType == EElasticity3D)
        {
            nStateVariables = 6;
        }
        TPZCompMesh *cmesh = fDistributedFlux.operator->();
        for(auto it : fMaterialIds)
        {
            TPZNullMaterial *nullmat = new TPZNullMaterial(it);
            nullmat->SetDimension(dim);
            nullmat->SetNStateVariables(nStateVariables);
            cmesh->InsertMaterialObject(nullmat);
        }
        cmesh = fAverageSolution.operator->();
        for(auto it : fMaterialIds)
        {
            TPZNullMaterial *nullmat = new TPZNullMaterial(it);
            nullmat->SetDimension(dim);
            nullmat->SetNStateVariables(nStateVariables);
            cmesh->InsertMaterialObject(nullmat);
        }
        cmesh = fCMeshLagrange.operator->();
        for(auto it : fMaterialIds)
        {
            TPZNullMaterial *nullmat = new TPZNullMaterial(it);
            nullmat->SetDimension(dim);
            nullmat->SetNStateVariables(nStateVariables);
            cmesh->InsertMaterialObject(nullmat);
        }
        cmesh = fCMeshConstantPressure.operator->();
        for(auto it : fMaterialIds)
        {
            TPZNullMaterial *nullmat = new TPZNullMaterial(it);
            nullmat->SetDimension(dim);
            nullmat->SetNStateVariables(nStateVariables);
            cmesh->InsertMaterialObject(nullmat);
        }
    }

}



/// Create all data structures for the computational mesh
void TPZMHMixedMeshControl::BuildComputationalMesh(bool usersubstructure)
{
    if (fpOrderInternal == 0 || fpOrderSkeleton == 0) {
        DebugStop();
    }
    InsertPeriferalMaterialObjects();
    CreateHDivMHMMesh();
    
#ifdef PZDEBUG
    if (fFluxMesh->Dimension() != fGMesh->Dimension()) {
        DebugStop();
    }
#endif
    
    CreatePressureMHMMesh();
    if(fNState > 1)
    {
        CreateRotationMesh();
    }
    
    
    CreateAverageSolutionMeshes();
    
    fLagrangeAveragePressure = true;
    CreateLagrangeMultiplierMesh();

    SetLagrangeMultiplierLevels();

    CreateHDivPressureMHMMesh();
    std::cout << "Total number of equations " << fCMesh->Solution().Rows() << std::endl;
    fGlobalSystemWithLocalCondensationSize = fCMesh->NEquations();
    fGlobalSystemSize = fCMesh->Solution().Rows();

    
#ifdef PZDEBUG
    if(1)
    {
        std::ofstream out1("HDivMesh.txt");
        fFluxMesh->Print(out1);
        std::ofstream out2("PressureMesh.txt");
        fPressureFineMesh->Print(out2);
        std::ofstream out3("DistributedFluxElementMesh.txt");
        fDistributedFlux->Print(out3);
        std::ofstream out4("AverageSolutionElementMesh.txt");
        fAverageSolution->Print(out4);
        std::ofstream out5("DistributedFluxMHMesh.txt");
        fCMeshLagrange->Print(out5);
        std::ofstream out6("AverageSolutionMHMesh.txt");
        fCMeshConstantPressure->Print(out6);
        std::ofstream out7("MultiPhysicsBeforeCondensing.txt");
        fCMesh->Print(out7);

    }
#endif
    
    fCMesh->ComputeNodElCon();
#ifdef PZDEBUG
    CheckMeshConsistency();
#endif
    
    if (usersubstructure) {
        HideTheElements();
    }
    fNumeq = fCMesh->NEquations();
    
#ifdef PZDEBUG
    {
        std::ofstream out("CondensedCMesh.txt");
        fCMesh->Print(out);
    }
#endif
#ifdef PZDEBUG2
    {
        int64_t nel = fCMesh->NElements();
        for(int64_t el = 0; el<nel; el++)
        {
            TPZCompEl *cel = fCMesh->Element(el);
            TPZSubCompMesh *sub = dynamic_cast<TPZSubCompMesh *>(cel);
            if(sub)
            {
                std::stringstream sout;
                sout << "submesh_" << el << ".vtk";
                std::ofstream file(sout.str());
                TPZVTKGeoMesh::PrintCMeshVTK(sub, file,true);
            }
        }

    }
#endif
}


void TPZMHMixedMeshControl::CreateHDivMHMMesh()
{
    TPZCompMesh * cmeshHDiv = fFluxMesh.operator->();
    cmeshHDiv->SetName("FluxMesh");
    InsertPeriferalHdivMaterialObjects();

    AdjustBoundaryElements();

    CreateInternalFluxElements();
    
#ifdef PZDEBUG
    if(0)
    {
        std::ofstream out("FluxmeshBeforeSkeleton.txt");
        fFluxMesh->Print(out);
    }
#endif
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        sout << "Subdomain indices " << fGeoToMHMDomain << std::endl;
        sout << "Recognized domains ";
        for (auto it = fMHMtoSubCMesh.begin(); it != fMHMtoSubCMesh.end(); it++) {
            sout << it->first << ' ';
        }
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    // Create the skeleton flux elements (boundary and internal)
    CreateSkeleton();
    
    // restrain the elements with matid fWrapMatId to the skeleton elements
    // will do nothing in hybrid hdiv meshes
    RestrainToSkeleton();
    
    if (fHdivmaismais) {
        int meshdim = fGMesh->Dimension();
        int64_t nel = cmeshHDiv->ElementVec().NElements();
        for (int64_t el = 0; el < nel; el++) {
            TPZCompEl *cel = cmeshHDiv->ElementVec()[el];
            TPZGeoEl *gel = cel->Reference();
            if (gel->Dimension() == meshdim) {
                int side = gel->NSides() - 1;
                TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cel);
                intel->SetSideOrder(side, fpOrderInternal + fHdivmaismais);//seta ordem +hdivmais
                intel->SetPreferredOrder(fpOrderInternal + fHdivmaismais);
            }
        }
        cmeshHDiv->ExpandSolution();
    }
    
    cmeshHDiv->ExpandSolution();
#ifdef PZDEBUG
    if(1)
    {
        fFluxMesh->ComputeNodElCon();
        std::ofstream outmesh("MixedMeshControl_HDivMesh.txt");
        Print(outmesh);
        std::ofstream outvtk("MixedMesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(fGMesh, outvtk,true);
    }
#endif
    return;
}

/// Insert the necessary H(div) material objects to create the flux mesh
void TPZMHMixedMeshControl::InsertPeriferalHdivMaterialObjects()
{
    TPZGeoMesh *gmesh = fGMesh.operator->();
    int meshdim = gmesh->Dimension();
    TPZCompMesh * cmeshHDiv = fFluxMesh.operator->();
    cmeshHDiv->SetDimModel(meshdim);
    if(fMaterialIds.size() == 0) DebugStop();
    TPZMaterial *matl2;
    for (auto item:fMaterialIds) {
        TPZNullMaterial *mat = new TPZNullMaterial(item);
        matl2 = mat;
        mat->SetNStateVariables(fNState);
        mat->SetDimension(meshdim);
        cmeshHDiv->InsertMaterialObject(matl2);
    }
    for (auto item:fMaterialBCIds) {
        TPZFNMatrix<1,STATE> val1(fNState,fNState,0.),val2(fNState,1,0.);
        TPZBndCond *bc = matl2->CreateBC(matl2, item, 0, val1, val2);
        cmeshHDiv->InsertMaterialObject(bc);
    }
    {
        TPZNullMaterial *nullmat = new TPZNullMaterial(fSkeletonMatId);
        nullmat->SetDimension(meshdim-1);
        nullmat->SetNStateVariables(fNState);
        cmeshHDiv->InsertMaterialObject(nullmat);
        if(fHybridizeSkeleton)
        {
            TPZNullMaterial *nullmat = new TPZNullMaterial(fSecondSkeletonMatId);
            nullmat->SetDimension(meshdim-1);
            nullmat->SetNStateVariables(fNState);
            cmeshHDiv->InsertMaterialObject(nullmat);
        }
    }
    {
        TPZNullMaterial *nullmat = new TPZNullMaterial(fWrapMatId);
        nullmat->SetDimension(meshdim-1);
        nullmat->SetNStateVariables(fNState);
        cmeshHDiv->InsertMaterialObject(nullmat);
    }
}



void TPZMHMixedMeshControl::CreatePressureMHMMesh()
{
    TPZGeoMesh * gmesh = fGMesh.operator->();
    gmesh->ResetReference();

    // the pressure mesh should be empty when calling this method
    int64_t npressureconnects = fPressureFineMesh->NConnects();
    if(npressureconnects != 0){
        DebugStop();
    }

    int porder = fpOrderInternal;
    // create and organize the pressure mesh
    // the pressure mesh is composed of discontinuous H1 elements
    TPZCompMesh * cmeshPressure = fPressureFineMesh.operator->();
    gmesh->ResetReference();
    cmeshPressure->SetName("PressureMesh");
    cmeshPressure->SetDimModel(gmesh->Dimension());
    if(porder+fHdivmaismais > 0)
    {
        cmeshPressure->SetAllCreateFunctionsContinuous();
    }
    else
    {
        cmeshPressure->SetAllCreateFunctionsDiscontinuous(); //AQUI
    }
    cmeshPressure->ApproxSpace().CreateDisconnectedElements(true);
    cmeshPressure->SetDefaultOrder(porder + fHdivmaismais);
    int meshdim = cmeshPressure->Dimension();
    // generate elements for all material ids of meshdim
    std::set<int> matids;
    for (auto it:fMaterialIds) {
        TPZMaterial *mat = cmeshPressure->FindMaterial(it);
        if (mat && mat->Dimension() == meshdim) {
            matids.insert(it);
        }
    }
    cmeshPressure->AutoBuild(matids);
    fPressureFineMesh->ExpandSolution();
    
    if(0)
    {
        std::ofstream out("PressureFineMesh.txt");
        fPressureFineMesh->Print(out);
    }

    
    
    // the lagrange multiplier level is set to one
    int64_t nc = cmeshPressure->NConnects();
    for (int64_t ic=0; ic<nc; ic++) {
        cmeshPressure->ConnectVec()[ic].SetLagrangeMultiplier(1);
    }
    // associate the connects with the proper subdomain
    gmesh->ResetReference();
    int64_t nel = cmeshPressure->NElements();
    for (int64_t el=0; el<nel; el++)
    {
        TPZCompEl *cel = cmeshPressure->Element(el);
#ifdef PZDEBUG
        if (! cel) {
            DebugStop();
        }
#endif
        // if a computational element was created outside the range of material ids
        // something very strange happened...
        TPZGeoEl *gel = cel->Reference();
        if(fMaterialIds.find (gel->MaterialId()) == fMaterialIds.end())
        {
            DebugStop();
        }
        int domain = fGeoToMHMDomain[gel->Index()];
#ifdef PZDEBUG
        if (domain == -1) {
            DebugStop();
        }
#endif

        SetSubdomain(cel, domain);
    }
    
    if(1){
        std::ofstream out("PressureFineMesh2.txt");
        fPressureFineMesh->Print(out);
    }
    
    return;
}


void TPZMHMixedMeshControl::CreateRotationMesh()
{
    TPZGeoMesh * gmesh = fGMesh.operator->();
    gmesh->ResetReference();
    int porder = fpOrderInternal;
    TPZCompMesh * cmeshRotation = fRotationMesh.operator->();
    gmesh->ResetReference();
    cmeshRotation->SetName("RotationMesh");
    cmeshRotation->SetDimModel(gmesh->Dimension());
    cmeshRotation->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
    cmeshRotation->SetDefaultOrder(porder);
    
    int meshdim = cmeshRotation->Dimension();
    std::set<int> matids;
    for (auto it:fMaterialIds) {
        TPZMaterial *mat = cmeshRotation->FindMaterial(it);
        if (mat && mat->Dimension() == meshdim) {
            matids.insert(it);
        }
    }
    cmeshRotation->AutoBuild(matids);
    fRotationMesh->ExpandSolution();
    
    if(0)
    {
        std::ofstream out("RotationsMesh.txt");
        fPressureFineMesh->Print(out);
    }

    
    int64_t nel = cmeshRotation->NElements();
    for(int64_t i=0; i<nel; i++){
        TPZCompEl *cel = cmeshRotation->ElementVec()[i];
        TPZCompElDisc *celdisc = dynamic_cast<TPZCompElDisc *>(cel);
        if(!celdisc) DebugStop();
        if(celdisc && celdisc->Reference()->Dimension() != meshdim)
        {
            DebugStop();
        }
        celdisc->SetTotalOrderShape();
        celdisc->SetFalseUseQsiEta();
    }
    
    int64_t nc = cmeshRotation->NConnects();
    for (int64_t ic=0; ic<nc; ic++) {
        cmeshRotation->ConnectVec()[ic].SetLagrangeMultiplier(1);
    }
    gmesh->ResetReference();
    
    for (int64_t el=0; el<nel; el++)
    {
        TPZCompEl *cel = cmeshRotation->Element(el);
#ifdef PZDEBUG
        if (! cel) {
            DebugStop();
        }
#endif
        TPZGeoEl *gel = cel->Reference();
        if(fMaterialIds.find (gel->MaterialId()) == fMaterialIds.end())
        {
            continue;
        }
#ifdef PZDEBUG
        if (fGeoToMHMDomain[gel->Index()] == -1) {
            DebugStop();
        }
#endif
        
        SetSubdomain(cel, fGeoToMHMDomain[gel->Index()]);
    }
    
    
    return;
}

/// Insert the necessary Pressure material objects to create the flux mesh
void TPZMHMixedMeshControl::InsertPeriferalPressureMaterialObjects()
{
    TPZCompMesh * cmeshPressure = fPressureFineMesh.operator->();
    
    
    for (auto it = fMaterialIds.begin(); it != fMaterialIds.end(); it++)
    {
        int matid = *it;
        if (cmeshPressure->MaterialVec().find(matid) == cmeshPressure->MaterialVec().end())
        {
            TPZNullMaterial *matl2 = new TPZNullMaterial((matid));
            matl2->SetNStateVariables(fNState);
            matl2->SetDimension(fGMesh->Dimension());
            cmeshPressure->InsertMaterialObject(matl2);
        }
        else
        {
            DebugStop();
        }
    }
    if(fPressureSkeletonMatId != 0)
    {
        if (fPressureFineMesh->FindMaterial(fPressureSkeletonMatId) != 0) {
            DebugStop();
        }
        TPZNullMaterial *mathybrid = new TPZNullMaterial(fPressureSkeletonMatId);
        mathybrid->SetDimension(fGMesh->Dimension()-1);
        mathybrid->SetNStateVariables(fNState);
        //fPressureFineMesh->InsertMaterialObject(mathybrid);
    }
    else
    {
        DebugStop();
    }

}

/// Insert the necessary Rotation material objects to create the flux mesh
void TPZMHMixedMeshControl::InsertPeriferalRotationMaterialObjects()
{
    TPZCompMesh * cmeshRotation = fRotationMesh.operator->();
    
    for (auto it = fMaterialIds.begin(); it != fMaterialIds.end(); it++)
    {
        int matid = *it;
        if (cmeshRotation->MaterialVec().find(matid) == cmeshRotation->MaterialVec().end())
        {
            TPZNullMaterial *matl2 = new TPZNullMaterial((matid));
            if(fNState == 2)
            {
                matl2->SetNStateVariables(1);
            } else if(fNState == 3)
            {
                matl2->SetNStateVariables(fNState);
            }
            matl2->SetDimension(fGMesh->Dimension());
            cmeshRotation->InsertMaterialObject(matl2);
        }
        else
        {
            DebugStop();
        }
    }
    
}

/// Put the pointers to the meshes in a vector
void TPZMHMixedMeshControl::GetMeshVec(TPZVec<TPZCompMesh *> &cmeshes)
{
    int nmeshes = 6;
    if(fNState > 1) nmeshes++;
    cmeshes.resize(nmeshes);
    cmeshes[0] = fFluxMesh.operator->();
    cmeshes[1] = fPressureFineMesh.operator->();
    if(fNState == 1)
    {
        cmeshes[2] = fDistributedFlux.operator->();
        cmeshes[3] = fAverageSolution.operator->();
        cmeshes[4] = this->fCMeshLagrange.operator->();
        cmeshes[5] = this->fCMeshConstantPressure.operator->();
    }
    else{
        cmeshes[2] = fRotationMesh.operator->();
        cmeshes[3] = fDistributedFlux.operator->();
        cmeshes[4] = fAverageSolution.operator->();
        cmeshes[5] = this->fCMeshLagrange.operator->();
        cmeshes[6] = this->fCMeshConstantPressure.operator->();
    }

}

/// Put the pointers to the meshes in a vector
void TPZMHMixedMeshControl::GetMeshVec(TPZVec<TPZAutoPointer<TPZCompMesh >> &cmeshes)
{
    int nmeshes = 6;
    if(fNState > 1) nmeshes++;
    cmeshes.resize(nmeshes);
    cmeshes[0] = fFluxMesh;
    cmeshes[1] = fPressureFineMesh;
    if(fNState == 1)
    {
        cmeshes[2] = fDistributedFlux;
        cmeshes[3] = fAverageSolution;
        cmeshes[4] = this->fCMeshLagrange;
        cmeshes[5] = this->fCMeshConstantPressure;
    }
    else{
        cmeshes[2] = fRotationMesh;
        cmeshes[3] = fDistributedFlux;
        cmeshes[4] = fAverageSolution;
        cmeshes[5] = this->fCMeshLagrange;
        cmeshes[6] = this->fCMeshConstantPressure;
    }

}

/// return the pointers to the meshes in a vector
TPZVec<TPZAutoPointer<TPZCompMesh> > TPZMHMixedMeshControl::GetMeshes()
{
    TPZVec<TPZAutoPointer<TPZCompMesh>> result;
    GetMeshVec(result);
    return result;
}



void TPZMHMixedMeshControl::CreateHDivPressureMHMMesh()
{
    TPZManVector<TPZCompMesh * > cmeshes;
    GetMeshVec(cmeshes);

    if(0)
    {
        std::ofstream out("PressureMesh_MultiPhis.txt");
        cmeshes[1] ->Print(out);
        
        std::ofstream out2("FluxMesh_MultiPhis.txt");
        cmeshes[0] ->Print(out2);
    }
    
    
    TPZGeoMesh *gmesh = fGMesh.operator->();
    if(!gmesh)
    {
        std::cout<< "Geometric mesh doesn't exist" << std::endl;
        DebugStop();
    }
    int dim = gmesh->Dimension();
    gmesh->ResetReference();
    // Multiphysics mesh
    TPZCompMesh * MixedFluxPressureCmesh = fCMesh.operator->();
    MixedFluxPressureCmesh->SetDimModel(dim);
    MixedFluxPressureCmesh->SetAllCreateFunctionsMultiphysicElem();
    
    TPZVec<int64_t> gelindices;
    this->GetGeometricElementPartition(gelindices);
    
    fCMesh->BuildMultiphysicsSpace(cmeshes,gelindices);

#ifdef PZDEBUG
    if(0)
    {
        std::ofstream out2("gmesh.txt");
        gmesh->Print(out2);
        std::ofstream out3("HDivMesh.txt");
        fFluxMesh->Print(out3);
        std::ofstream out4("PressureMesh.txt");
        fPressureFineMesh->Print(out4);
    }
#endif
    

    // populate the connect to subdomain data structure for the multiphysics mesh
    JoinSubdomains(cmeshes, MixedFluxPressureCmesh);
    
    // Transferindo para a multifisica
//    TPZBuildMultiphysicsMesh::AddElements(meshvector, MixedFluxPressureCmesh);
//    TPZBuildMultiphysicsMesh::AddConnects(meshvector, MixedFluxPressureCmesh);
    // copy the solution of the atomic meshes to the multiphysics mesh
    TPZBuildMultiphysicsMesh::TransferFromMeshes(cmeshes, MixedFluxPressureCmesh);

    if(fHybridizeSkeleton) DebugStop();
#ifdef PZDEBUG
    if(0)
    {
        MixedFluxPressureCmesh->ComputeNodElCon();
        std::ofstream file("cmeshmphys.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(MixedFluxPressureCmesh, file,true);
        std::ofstream out("cmeshmphys.txt");
        MixedFluxPressureCmesh->Print(out);
    }
#endif
    MixedFluxPressureCmesh->CleanUpUnconnectedNodes();
    
#ifdef PZDEBUG
    if(0)
    {
        std::ofstream out("multiphysics.txt");
        MixedFluxPressureCmesh->Print(out);
    }
#endif
    
    return;
    
}

void TPZMHMixedMeshControl::HideTheElements()
{
    char ExternalLagrangeLevel = 11;
    if (!fLagrangeAveragePressure) {
        ExternalLagrangeLevel = 9;
    }
    typedef std::set<int64_t> TCompIndexes;
    std::map<int64_t, TCompIndexes> ElementGroups;
    TPZGeoMesh *gmesh = fCMesh->Reference();
    gmesh->ResetReference();
    int dim = gmesh->Dimension();
    fCMesh->LoadReferences();
    int64_t nel = fCMesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = fCMesh->Element(el);
        int64_t domain = WhichSubdomain(cel);
        if (domain == -1) {
            continue;
        }
        ElementGroups[domain].insert(el);
    }
    if (ElementGroups.size() <= 10)
    {
        std::cout << "Number of element groups " << ElementGroups.size() << std::endl;
        std::map<int64_t,TCompIndexes>::iterator it;
        for (it=ElementGroups.begin(); it != ElementGroups.end(); it++) {
            std::cout << "Group " << it->first << " group size " << it->second.size() << std::endl;
            std::cout << " elements ";
            std::set<int64_t>::iterator its;
            for (its = it->second.begin(); its != it->second.end(); its++) {
                std::cout << *its << "|" << fCMesh->Element(*its)->Reference()->Index() << " ";
            }
            std::cout << std::endl;
        }
    }
    
    std::map<int64_t,int64_t> submeshindices;
    TPZCompMeshTools::PutinSubmeshes(fCMesh.operator->(), ElementGroups, submeshindices, ExternalLagrangeLevel);
    std::cout << "After putting in substructures\n";
    fMHMtoSubCMesh = submeshindices;
    fCMesh->ComputeNodElCon();
    fCMesh->CleanUpUnconnectedNodes();
    
    GroupandCondenseElements();

    std::cout << "Finished substructuring\n";
}

void TPZMHMixedMeshControl::HideTheElements(int KeepOneLagrangian)
{
    typedef std::set<int64_t> TCompIndexes;
    std::map<int64_t, TCompIndexes> ElementGroups;
    TPZGeoMesh *gmesh = fCMesh->Reference();
    gmesh->ResetReference();
    int dim = gmesh->Dimension();
    fCMesh->LoadReferences();
    int64_t nel = fCMesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = fCMesh->Element(el);
        int64_t domain = WhichSubdomain(cel);
        if (domain == -1) {
            continue;
        }
        ElementGroups[domain].insert(el);
    }
    if (ElementGroups.size() <= 10)
    {
        std::cout << "Number of element groups " << ElementGroups.size() << std::endl;
        std::map<int64_t,TCompIndexes>::iterator it;
        for (it=ElementGroups.begin(); it != ElementGroups.end(); it++) {
            std::cout << "Group " << it->first << " group size " << it->second.size() << std::endl;
            std::cout << " elements ";
            std::set<int64_t>::iterator its;
            for (its = it->second.begin(); its != it->second.end(); its++) {
                std::cout << *its << "|" << fCMesh->Element(*its)->Reference()->Index() << " ";
            }
            std::cout << std::endl;
        }
    }
    
    std::map<int64_t,int64_t> submeshindices;
    for (auto elgroup : ElementGroups) {
        int64_t submeshindex;
        TPZCompMeshTools::PutinSubmeshes(fCMesh.operator->(), elgroup.second, submeshindex, KeepOneLagrangian);
        submeshindices[elgroup.first] = submeshindex;
    }
//    TPZCompMeshTools::PutinSubmeshes(fCMesh.operator->(), ElementGroups, submeshindices, KeepOneLagrangian);
    std::cout << "After putting in substructures\n";
    fMHMtoSubCMesh = submeshindices;
    fCMesh->ComputeNodElCon();
    fCMesh->CleanUpUnconnectedNodes();
    
    GroupandCondenseElements();

    std::cout << "Finished substructuring\n";
}


/// print the data structure
void TPZMHMixedMeshControl::Print(std::ostream &out)
{
    
    /// geometric mesh used to create the computational mesh
    if (fGMesh)
    {
        out << "******************* GEOMETRIC MESH *****************\n";
        fGMesh->Print(out);
    }
    
    /// computational mesh to contain the pressure elements
    // this mesh is the same as fCMesh if there are no lagrange multipliers assocated with the average pressure
    if (fFluxMesh)
    {
        out << "******************* FLUX MESH *****************\n";
        fFluxMesh->Print(out);
    }
    
    /// computational mesh to contain the pressure elements
    // this mesh is the same as fCMesh if there are no lagrange multipliers assocated with the average pressure
    if (fPressureFineMesh)
    {
        out << "******************* PRESSURE MESH *****************\n";
        fPressureFineMesh->Print(out);
    }
    
    
    /// computational MHM mesh being built by this class
    if (fCMesh)
    {
        out << "******************* COMPUTATIONAL MESH *****************\n";
        fCMesh->Print(out);
    }
    
    /// computational mesh to represent the constant states
    if (fCMeshLagrange)
    {
        out << "******************* LAGRANGE MULTIPLIER MESH *****************\n";
        fCMeshLagrange->Print(out);
    }
    
    /// computational mesh to represent the constant states
    if (fCMeshConstantPressure)
    {
        out << "******************* CONSTANTE PRESSURE MESH *****************\n";
        fCMeshConstantPressure->Print(out);
    }
    
    /// material id associated with the skeleton elements
    out << "Skeleton Mat Id " <<  fSkeletonMatId << std::endl;
    
    /// material id associated with the lagrange multiplier elements
    out << "Lagrange mat id left - right " <<  fLagrangeMatIdLeft << " - " << fLagrangeMatIdRight << std::endl;
    
    /// interpolation order of the internal elements
    out << "Internal polynomial order " <<  fpOrderInternal << std::endl;
    
    /// interpolation order of the skeleton elements
    out << "Skeleton polynomial order " << fpOrderSkeleton << std::endl;
    
    /// indices of the geometric elements which define the skeleton mesh
    {
        out << "Geometric element indices of the coarse mesh ";
        for (std::map<int64_t,int64_t>::iterator it= fMHMtoSubCMesh.begin(); it != fMHMtoSubCMesh.end(); it++) {
            out << it->first << " " << it->second << " ";
        }
        out << std::endl;
    }
    /// indices of the skeleton elements and their left/right elements of the skeleton mesh
    out << "Skeleton element indices with associated left and right coarse element indices\n";
    {
        std::map<int64_t, std::pair<int64_t,int64_t> >::iterator it;
        for (it = fInterfaces.begin(); it != fInterfaces.end(); it++) {
            out << "skel index " << it->first << " Left Right indices " << it->second.first << " " << it->second.second << std::endl;
        }
    }
    
    /// flag to determine whether a lagrange multiplier is included to force zero average pressures in the subdomains
    /**
     * when imposing average pressure to be zero, a multiphysics mesh is created
     */
    out << "Will generate a constant pressure mesh " <<  fLagrangeAveragePressure << std::endl;
    /// subdonain indices of the connects
    out << "Subdomain indices of the connects by mesh\n";
    for (auto it = fConnectToSubDomainIdentifier.begin(); it != fConnectToSubDomainIdentifier.end(); it++) {
        TPZCompMesh *cmesh = it->first;
        TPZManVector<int64_t> &cvec = it->second;
        out << "Computational mesh pointer " << (void *) cmesh << std::endl;
        for (int64_t i=0; i<cvec.size(); i++) {
            if (i && !(i%20)) {
                out << std::endl;
            }
            out << "(" << i << "->" << cvec[i] << ") ";
        }
        out << std::endl;
    }
    
}

static void SetSubMatid(TPZGeoEl *gel, int matid)
{
    int nsub = gel->NSubElements();
    for (int isub=0; isub < nsub; isub++) {
        TPZGeoEl *subel = gel->SubElement(isub);
        subel->SetMaterialId(matid);
        SetSubMatid(subel, matid);
    }
}
static void SetAllMatid(TPZGeoEl *gel, int matid)
{
    gel->SetMaterialId(matid);
    TPZGeoEl *fat = gel->Father();
    while (fat) {
        fat->SetMaterialId(matid);
        fat = fat->Father();
    }
    SetSubMatid(gel, matid);
}

void TPZMHMixedMeshControl::HybridizeSkeleton(int skeletonmatid, int pressurematid)
{

    fGMesh->ResetReference();
    int64_t nskel = fInterfaces.size();
    std::map<int64_t, std::pair<int64_t, int64_t> >::iterator it;
    TPZStack<TPZInterpolatedElement*> fluxorig;
    TPZStack<TPZInterpolatedElement *> fluxsecond;
    TPZStack<TPZInterpolatedElement *> pressure;

    fFluxMesh->LoadReferences();
    // build the fluxorig datastructure : contains the original flux elements
    // loop over the skeleton elements
    for (it = fInterfaces.begin(); it != fInterfaces.end(); it++) {
        TPZGeoEl *gel = fGMesh->Element(it->first);
        // skip the boundary elements
        if (it->first == it->second.second) {
            continue;
        }
        if (gel->MaterialId() != skeletonmatid) {
            continue;
        }
        int side = gel->NSides()-1;
        TPZInterpolatedElement *orig = dynamic_cast<TPZInterpolatedElement *>(gel->Reference());
        fluxorig.Push(orig);
    }
    fGMesh->ResetReference();
    
    fPressureFineMesh->SetDefaultOrder(fpOrderSkeleton);
    // first create a second flux element and a pressure element on top of the existing skeleton element
    // loop over the skeleton elements
    for (it = fInterfaces.begin(); it != fInterfaces.end(); it++) {
        TPZGeoEl *gel = fGMesh->Element(it->first);
        // skip the boundary elements
        if (it->first == it->second.second) {
            continue;
        }
        if (gel->MaterialId() != skeletonmatid) {
            continue;
        }
        int side = gel->NSides()-1;
        TPZGeoElSide gelside(gel,side);
        TPZGeoElBC skeleton2(gelside,fSecondSkeletonMatId);
        fFluxMesh->SetAllCreateFunctionsHDiv();
        // create a flux boundary element
        int64_t indexflux;
        fFluxMesh->CreateCompEl(skeleton2.CreatedElement(), indexflux);
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(fFluxMesh->Element(indexflux));
        SetSubdomain(intel, -1);
#ifdef PZDEBUG
        if (intel->NConnects() != 1) {
            DebugStop();
        }
#endif
        // the side orientation of the boundary fluxes is +1 - we will need to fix the elements as well!
        if (it->second.first < it->second.second) {
            // totototo
            intel->SetSideOrient(side, 1);
        }
        else
        {
            intel->SetSideOrient(side, 1);
        }
        skeleton2.CreatedElement()->ResetReference();
        
        // create a dim-1 dimensional pressure element
        TPZGeoElBC pressuregel(gelside,pressurematid);
        // this will be changed to the pressure mesh
//        fFluxMesh->SetAllCreateFunctionsContinuous();
        fPressureFineMesh->SetAllCreateFunctionsContinuous();
        int64_t indexpressure;
//        fFluxMesh->CreateCompEl(pressuregel.CreatedElement(), indexpressure);
        fPressureFineMesh->CreateCompEl(pressuregel.CreatedElement(), indexpressure);
        TPZInterpolatedElement *presel = dynamic_cast<TPZInterpolatedElement *>(fPressureFineMesh->Element(indexpressure));
        SetSubdomain(presel, -1);
        // set the lagrange multiplier to the highest level
        int nc = presel->NConnects();
        for (int ic=0; ic<nc; ic++) {
            TPZConnect &c = presel->Connect(ic);
            c.SetLagrangeMultiplier(3);
        }
        // This can only be done after all flux connects have been created!!!
//        SetSubdomain(presel, -1);
        pressuregel.CreatedElement()->ResetReference();
        fluxsecond.Push(intel);
        pressure.Push(presel);
    }
    fFluxMesh->LoadReferences();
    
    // switch the connect dependency around AND fix the side orientation
    int64_t count = 0;
    for (it = fInterfaces.begin(); it != fInterfaces.end(); it++) {
        TPZGeoEl *gel = fGMesh->Element(it->first);
        // skip the boundary elements
        if (it->first == it->second.second) {
            continue;
        }
        if (gel->MaterialId() != skeletonmatid) {
            continue;
        }
        TPZCompEl *cel = gel->Reference();
        if (!cel || cel != fluxorig[count]) {
            DebugStop();
        }
        std::map<int64_t,std::list<TPZCompElSide> > connectedmap;
        // connected elements will compute all the subelements of the skeleton element to left and right
        ConnectedElements(it->first, it->second, connectedmap);
        if (connectedmap.size() != 2) {
            DebugStop();
        }
        // identify the domains left and right of the skeleton element
        int64_t origdepindex = fluxorig[count]->ConnectIndex(0);
        int64_t newdepindex = fluxsecond[count]->ConnectIndex(0);
        std::list<TPZCompElSide> &updatelist = connectedmap[it->second.second];
        for (std::list<TPZCompElSide>::iterator itlist = updatelist.begin(); itlist != updatelist.end(); itlist++)
        {
            TPZCompElSide smallCompElSide = *itlist;
            TPZInterpolatedElement *smallel = dynamic_cast<TPZInterpolatedElement *>(smallCompElSide.Element());
            if (smallel->NSideConnects(smallCompElSide.Side()) != 1) {
                DebugStop();
            }
            TPZConnect &c = smallel->SideConnect(0, smallCompElSide.Side());
            TPZConnect::TPZDepend *dep = c.FirstDepend();
            if(dep->fDepConnectIndex != origdepindex)
            {
                DebugStop();
            }
            dep->fDepConnectIndex = newdepindex;
            // Set the side orientation to +1
            smallel->SetSideOrient(smallCompElSide.Side(), 1);
        }
        SetSubdomain(fluxorig[count],it->second.first);
        SetSubdomain(fluxsecond[count], it->second.second);
        count++;
    }
    fGMesh->ResetReference();
    // switch the reference index of the original flux element and the pressure element
    int64_t numflux = fluxorig.size();
    for (int iflux=0; iflux<numflux; iflux++) {
        int64_t fluxgelindex = fluxorig[iflux]->Reference()->Index();
        int fluxmat = fluxorig[iflux]->Reference()->MaterialId();
        int pressmat = pressure[iflux]->Reference()->MaterialId();
        int64_t pressuregelindex = pressure[iflux]->Reference()->Index();
        fluxorig[iflux]->SetReference(pressuregelindex);
        pressure[iflux]->SetReference(fluxgelindex);
        SetAllMatid(fluxorig[iflux]->Reference(), fluxmat);
        SetAllMatid(pressure[iflux]->Reference(), pressmat);
    }

    // The connects of the pressure mesh are not to be condensed
    fFluxMesh->ExpandSolution();
    fPressureFineMesh->ExpandSolution();

}

// create the elements domain per domain with approximation spaces disconnected from each other
void TPZMHMixedMeshControl::CreateInternalFluxElements() {
    TPZGeoMesh *gmesh = fGMesh.operator->();
    int meshdim = gmesh->Dimension();
    TPZCompMesh *cmeshHDiv = fFluxMesh.operator->();
    gmesh->ResetReference();
    cmeshHDiv->LoadReferences();
    cmeshHDiv->SetDimModel(meshdim);
    cmeshHDiv->ApproxSpace().SetAllCreateFunctionsHDiv(meshdim);
    cmeshHDiv->SetDefaultOrder(fpOrderInternal);

    // Create MHM computational elements
    fConnectToSubDomainIdentifier[cmeshHDiv].Expand(10000);
    int64_t nelgeo = fGMesh->NElements();

    if (fGeoToMHMDomain.size() != nelgeo) DebugStop();

    TPZManVector<std::pair<int64_t, int64_t>> MHMOfEachGeoEl(nelgeo);
    for (int i = 0; i < nelgeo; i++) {
        MHMOfEachGeoEl[i] = std::make_pair(fGeoToMHMDomain[i], i);
    }

    std::sort(&MHMOfEachGeoEl[0], &MHMOfEachGeoEl[nelgeo - 1] + 1);

    int64_t previousMHMDomain = -1;
    int64_t firstElemInMHMDomain = -1;
    for (int i = 0; i < nelgeo; i++) {
        int64_t MHMDomain = MHMOfEachGeoEl[i].first;
        int64_t elIndex = MHMOfEachGeoEl[i].second;

        if (MHMDomain == -1) continue;

        if (MHMDomain != previousMHMDomain) {
            if (previousMHMDomain != -1) {
                // create a wrap mesh around the elements on the boundary of the
                // mhm domain
                for (int j = firstElemInMHMDomain; j < i; j++) {
                    fGMesh->Element(MHMOfEachGeoEl[j].second)->ResetReference();
                }
            }
            firstElemInMHMDomain = i;
            previousMHMDomain = MHMDomain;
        }

        // Create the flux element
        TPZGeoEl *gel = fGMesh->Element(elIndex);
        if (!gel || gel->HasSubElement()) continue;
        int64_t index;
        cmeshHDiv->CreateCompEl(gel, index);
        TPZCompEl *cel = cmeshHDiv->Element(index);
        // Associate the connects with the subdomain for the flux mesh
        SetSubdomain(cel, MHMDomain);
    }
    // Resets references of last MHM domain
    for (int j = firstElemInMHMDomain; j < nelgeo; j++) {
        fGMesh->Element(MHMOfEachGeoEl[j].second)->ResetReference();
    }
    fFluxMesh->ExpandSolution();
    fFluxMesh->ComputeNodElCon();
    int64_t nelflux = fFluxMesh->NElements();
    for (int64_t el = 0; el<nelflux; el++) {
        TPZCompEl *cel = fFluxMesh->Element(el);
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        TPZGeoEl *gel = cel->Reference();
        if(gel->Dimension() != meshdim) continue;
        int nfaces = gel->NSides(meshdim-1);
        for (int side = gel->NSides()-nfaces-1; side < gel->NSides()-1; side++) {
            TPZConnect &c = intel->SideConnect(0, side);
            if(c.HasDependency()) continue;
            if(c.NElConnected() == 1)
            {
                TPZGeoElBC gelbc(gel,side,fWrapMatId);
                int subdomain = WhichSubdomain(gel);
                cel->LoadElementReference();
                int64_t index;
                TPZCompEl *wrap = fFluxMesh->CreateCompEl(gelbc.CreatedElement(), index);
                SetSubdomain(wrap, subdomain);
                gel->ResetReference();
                gelbc.CreatedElement()->ResetReference();
            }
        }
    }
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        fFluxMesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
}

/// will create the elements on the skeleton
void TPZMHMixedMeshControl::CreateSkeletonOld()
{
    // comment this line or not to switch the type of skeleton elements
    int meshdim = fFluxMesh->Dimension();
    fFluxMesh->ApproxSpace().SetAllCreateFunctionsHDiv(meshdim);
    fGMesh->ResetReference();
    int order = fpOrderSkeleton;
    if (order <= 0) {
        DebugStop();
        order = 1;
    }
    // create the skeleton elements without applying the restraints of the elements of the subdomains
    fFluxMesh->SetDefaultOrder(order);
    std::map<int64_t, std::pair<int64_t, int64_t> >::iterator it = fInterfaces.begin();
    while (it != fInterfaces.end()) {
        int64_t elindex = it->first;

        if (elindex == it->second.second) {
            DebugStop();
            // boundary elements shouldn't be in fInterface map. they are removed from the map in method
            // TPZMHMixedMeshControl::AdjustBoundaryElements
        }
        TPZGeoEl *gel = fGMesh->ElementVec()[elindex];
        int64_t index;
        // create an element to model the flux between subdomains
        fFluxMesh->CreateCompEl(gel, index);
        TPZCompEl *cel = fFluxMesh->ElementVec()[index];
        int Side = gel->NSides()-1;
        TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
        SetSubdomain(cel, -1);

        if (it->second.first < it->second.second) {
          // set the flux orientation depending on the relative value of the
          // element ids
          intel->SetSideOrient(Side, 1);
        } else {
          intel->SetSideOrient(Side, -1);
        }
        // this element will not be put in a subdomain
        SetSubdomain(cel, -1);
        gel->ResetReference();

        it++;
    }
    // Apply restraints to the element/sides and the skeleton
    fFluxMesh->LoadReferences();
#ifdef PZDEBUG
    if(0)
    {
        std::ofstream out("FluxBeforeInterfaces.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(fFluxMesh.operator->(), out);
    }
#endif
#ifdef PZDEBUG
    // verify if any connect is restrained
    {
        int64_t nconnects = fFluxMesh->NConnects();
        for (int64_t ic = 0; ic<nconnects; ic++) {
            TPZConnect &c = fFluxMesh->ConnectVec()[ic];
            if (c.HasDependency()) {
                std::cout << "Connect index " << ic << " Has dependency\n";
            }
        }
    }
#endif
    it = fInterfaces.begin();
    while (it != fInterfaces.end()) {
        // the index of the map of the interface data structure is the index of the geometric element
        int64_t elindex = it->first;

        // data structure containing the list of connected elements by subdomain
        std::map<int64_t,std::list<TPZCompElSide> > subels;
        // find the connected computational elements in function of the index of the geometric element
        ConnectedElements(elindex, it->second, subels);

#ifdef LOG4CXX
        if(logger->isDebugEnabled())
        {
            std::stringstream sout;
            sout << "Interface elindex " << elindex << " left " << it->second.first << " right " << it->second.second << std::endl;
            sout << "Number of connected subdomains " << subels.size() << std::endl;
            for (std::map<int64_t,std::list<TPZCompElSide> >::iterator itlist = subels.begin(); itlist != subels.end(); itlist++)
            {
                sout << " SUBDOMAIN first " << itlist->first << "\n";
                std::list<TPZCompElSide> &lst = itlist->second;
                for (std::list<TPZCompElSide>::iterator its = lst.begin(); its != lst.end(); its++)
                {
                    sout << " GeoElSide " << its->Reference() << std::endl;
                }
            }
            LOGPZ_DEBUG(logger,sout.str())
        }
#endif
        TPZGeoEl *gel = fGMesh->ElementVec()[elindex];
        int side = gel->NSides()-1;
        TPZCompEl *cel = gel->Reference();
        // intel contains the flux skeleton element
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);

        for (std::map<int64_t,std::list<TPZCompElSide> >::iterator itlist = subels.begin(); itlist != subels.end(); itlist++) {
            std::list<TPZCompElSide> &lst = itlist->second;
            for (std::list<TPZCompElSide>::iterator its = lst.begin(); its != lst.end(); its++) {
                TPZInterpolatedElement *intelsub = dynamic_cast<TPZInterpolatedElement *>(its->Element());
                if (!intelsub) {
                    DebugStop();
                }
                intelsub->RestrainSide(its->Side(), intel, side);
            }
        }
        it++;
    }
    // this will renumber the connects that are now constrained
    fFluxMesh->ExpandSolution();
    fFluxMesh->SetDimModel(meshdim);
}



/// Create the interfaces between the pressure elements of dimension dim
void TPZMHMixedMeshControl::CreateMultiPhysicsInterfaceElements(int dim)
{
    std::pair<int,int> skelmatid(fSkeletonMatId,fSecondSkeletonMatId);
    CreateMultiPhysicsInterfaceElements(dim, fPressureSkeletonMatId, skelmatid);
}


/// Create the interfaces between the pressure elements of dimension dim
/// This method will create interface elements between the pressure element and
// neighbouring flux elements
void TPZMHMixedMeshControl::CreateMultiPhysicsInterfaceElements(int dim, int pressmatid, std::pair<int,int> hdivmatid)
{
    TPZCompMesh * cmeshPressure = fPressureFineMesh.operator->();
    
    // Computational multiplysics mesh
    TPZCompMesh * MixedFluxPressureCmesh = fCMesh.operator->();
    // load only dim dimensional elements
    MixedFluxPressureCmesh->Reference()->ResetReference();
    int64_t nel = MixedFluxPressureCmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = MixedFluxPressureCmesh->Element(el);
        if (!cel || !cel->Reference()) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (gel->Dimension() != dim) {
            continue;
        }
        cel->LoadElementReference();
    }

    // for each pressure element neighbour of a skeleton element create two interface elements
    nel = cmeshPressure->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = cmeshPressure->Element(el);
        if (!cel || !cel->Reference()) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (gel->Dimension() != dim || gel->MaterialId() != pressmatid) {
            continue;
        }
        if(!gel->Reference())
        {
            DebugStop();
        }
        TPZGeoElSide gelside(gel,gel->NSides()-1);
        // look for an fSkeletonWrapMatId element
        TPZGeoElSide neighbour = gelside.Neighbour();
        while (neighbour != gelside) {
            if (neighbour.Element()->MaterialId() == fWrapMatId) {
                break;
            }
            neighbour = neighbour.Neighbour();
        }
        // we ignore elements without neighbouring fSkeletonWrapMatId?
        if (neighbour == gelside) {
            continue;
        }
        TPZStack<TPZCompElSide> celstack;
        // This method will put all equallevel elements on the stack, including the reference to the current element
        // Reference to the current element will be put last on the stack
        neighbour = gelside.Neighbour();
        while (neighbour != gelside) {
            if (neighbour.Element()->MaterialId() == hdivmatid.first || neighbour.Element()->MaterialId() == hdivmatid.second) {
                celstack.Push(neighbour.Reference());
            }
            neighbour = neighbour.Neighbour();
        }
        if (celstack.size() == 1) {
            int matid = celstack[0].Element()->Reference()->MaterialId();
            if (matid == fSecondSkeletonMatId) {
                continue;
            }
        }
        if (celstack.size() != 2) {
            std::cout << "Error found only material ids ";
            for (int i=0; i<celstack.size(); i++) {
                std::cout << celstack[i].Element()->Reference()->MaterialId() << " ";
            }
            std::cout << std::endl;
            DebugStop();
        }
        TPZCompElSide celside = gelside.Reference();
        // create an interface between the multiphysics element associated with the pressure element
        // and the neighbouring elements
        TPZGeoElBC gbcleft(gelside, fLagrangeMatIdLeft);
        TPZGeoElBC gbcright(gelside, fLagrangeMatIdRight);
        int64_t index1,index2;
//        celside.Element()->Print();
//        celstack[0].Element()->Print();
//        celstack[1].Element()->Print();
        new TPZMultiphysicsInterfaceElement(*MixedFluxPressureCmesh,gbcleft.CreatedElement(),index1,celstack[0],celside);
        new TPZMultiphysicsInterfaceElement(*MixedFluxPressureCmesh,gbcright.CreatedElement(),index2,celstack[1],celside);
#ifdef LOG4CXX
        if(logger->isDebugEnabled())
        {
            std::stringstream sout;
            sout << "Created an interface skeleton element index " << index1 << " between mphys " << celstack[0].Element()->Index() << " and pressure mphyx " << celside.Element()->Index()
            << std::endl;
            sout << "Created an interface skeleton element index " << index2 << " between mphys " << celstack[1].Element()->Index() << " and pressure mphyx " << celside.Element()->Index();
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
    }
}

/// group and condense the elements
void TPZMHMixedMeshControl::GroupandCondenseElements()
{
    for (std::map<int64_t,int64_t>::iterator it=fMHMtoSubCMesh.begin(); it != fMHMtoSubCMesh.end(); it++) {
        TPZCompEl *cel = fCMesh->Element(it->second);
        TPZSubCompMesh *subcmesh = dynamic_cast<TPZSubCompMesh *>(cel);
        if (!subcmesh) {
            DebugStop();
        }
        TPZCompMeshTools::GroupElements(subcmesh);
        subcmesh->ComputeNodElCon();
        
#ifdef LOG4CXX2
        if(logger->isDebugEnabled())
        {
            std::stringstream sout;
            subcmesh->Print(sout);
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        // TODO: increment nelconnected of exterior connects
        // the average element pressure should not be condensed
        char keeplagrange = 6;
        bool keepmatrix = false;
        TPZCompMeshTools::CondenseElements(subcmesh, keeplagrange, keepmatrix);
        subcmesh->CleanUpUnconnectedNodes();
        int numthreads = 0;
        int preconditioned = 0;
#ifdef LOG4CXX2
        if(logger->isDebugEnabled())
        {
            std::stringstream sout;
            subcmesh->Print(sout);
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        TPZAutoPointer<TPZGuiInterface> guiInterface;
        subcmesh->SetAnalysisSkyline(numthreads, preconditioned, guiInterface);
        // @TODO change the lagrange multiplier levels in function of the equation numbers and then perform Saddle permute
    }
    
}

/// delete the pressure elements leaving the geometric mesh without pointing to the computational mesh
void TPZMHMixedMeshControl::DeletePressureElements()
{
    fGMesh->ResetReference();
    int64_t nel = fPressureFineMesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = fPressureFineMesh->Element(el);
        if (cel) {
            delete cel;
        }
    }
}

/// build the multi physics mesh (not at the finest geometric mesh level
void TPZMHMixedMeshControl::BuildMultiPhysicsMesh()
{
    if (fCMesh->NElements() != 0) {
        DebugStop();
    }
    fCMesh->SetAllCreateFunctionsMultiphysicElem();
    TPZMultiphysicsCompMesh *mphysics = dynamic_cast<TPZMultiphysicsCompMesh *>(fCMesh.operator->());
    int vecsize = 2;
    if(fNState > 1) vecsize = 3;
    TPZManVector<TPZCompMesh *> meshvec(vecsize);
    meshvec[0] = fFluxMesh.operator->();
    meshvec[1] = fPressureFineMesh.operator->();
    if(fNState > 1)
    {
        meshvec[2] = this->fRotationMesh.operator->();
    }
    TPZManVector<int64_t> shouldcreate(fGMesh->NElements(),0);
    std::set<int> matids;
    for (auto it : fCMesh->MaterialVec()) {
        matids.insert(it.first);
    }
    int64_t nel = fFluxMesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = fFluxMesh->Element(el);
        TPZGeoEl *gel = cel->Reference();
        // this means that all geometric elements associated with flux elements will generate a computational element
        if (matids.find(gel->MaterialId()) == matids.end()) {
            DebugStop();
        }
        if (cel) {
            shouldcreate[cel->Reference()->Index()] = 1;
        }
    }
    nel = fPressureFineMesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = fPressureFineMesh->Element(el);
        TPZGeoEl *gel = cel->Reference();
        if (matids.find(gel->MaterialId()) == matids.end()) {
            DebugStop();
        }
        if (cel) {
            shouldcreate[cel->Reference()->Index()] = 1;
        }
    }
    // define the intersection of the finest references
    nel = fGMesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = fGMesh->Element(el);
        if (shouldcreate[el])
        {
            TPZGeoEl *fat = gel->Father();
            while(fat)
            {
                if(shouldcreate[fat->Index()] == 1)
                {
                    shouldcreate[fat->Index()] = 0;
                }
                fat = fat->Father();
            }
        }
    }
    TPZStack<int64_t> gelindexes;
    for (int64_t el=0; el<nel; el++) {
        if (shouldcreate[el])
        {
            gelindexes.Push(el);
        }
    }
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << "Geometric indices for which we will create multiphysics elements" << std::endl;
        sout << gelindexes;
//        std::cout << sout.str() << std::endl;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    mphysics->BuildMultiphysicsSpace(meshvec,gelindexes);

}

/// create the distributed flux and average solution mesh
void TPZMHMixedMeshControl::CreateAverageSolutionMeshes()
{
    if(!fAverageSolution || !fDistributedFlux) DebugStop();
    int64_t nel_pressure = fPressureFineMesh->NElements();
    int meshdim = fGMesh->Dimension();
    TPZStack<int64_t> gelindexes;
    for (int el = 0; el<nel_pressure; el++) {
        TPZCompEl *cel = fPressureFineMesh->Element(el);
        TPZGeoEl *gel = cel->Reference();
        int dim = gel->Dimension();
        if(dim == meshdim) gelindexes.Push(gel->Index());
    }
    fGMesh->ResetReference();
    fAverageSolution->SetAllCreateFunctionsDiscontinuous();
    fAverageSolution->SetDefaultOrder(0);
    fAverageSolution->AutoBuild(gelindexes);
    TPZCompElDisc::SetTotalOrderShape(fAverageSolution.operator->());
    fGMesh->ResetReference();
    fDistributedFlux->SetAllCreateFunctionsDiscontinuous();
    fDistributedFlux->SetDefaultOrder(0);
    fDistributedFlux->AutoBuild(gelindexes);
    TPZCompElDisc::SetTotalOrderShape(fDistributedFlux.operator->());

    {
        int64_t nel = fAverageSolution->NElements();
        for (int64_t el=0; el<nel; el++) {
            TPZCompEl *cel = fAverageSolution->Element(el);
            TPZCompEl *cel2 = fDistributedFlux->Element(el);
            TPZGeoEl *gel = cel->Reference();
            int subdomain = fGeoToMHMDomain[gel->Index()];
            SetSubdomain(cel, subdomain);
            SetSubdomain(cel2, subdomain);
        }
    }

}

// restrain the connects to the skeleton and boundary elements
void TPZMHMixedMeshControl::RestrainToSkeleton()
{
    fGMesh->ResetReference();
    fFluxMesh->LoadReferences();
    
    std::set<int> targetmatid(fMaterialBCIds);
    targetmatid.insert(fSkeletonMatId);
    
    int64_t nelflux = fFluxMesh->NElements();
    TPZCompMesh *cmesh = fFluxMesh.operator->();
    for (int64_t el=0; el<nelflux; el++) {
        TPZCompEl *cel = fFluxMesh->Element(el);
        TPZGeoEl *gel = cel->Reference();
        int matid = gel->MaterialId();
        if(matid != fWrapMatId) continue;
        TPZGeoElSide gelside(gel);
        TPZGeoElSideAncestors ancestor(gelside);
        TPZGeoElSide neighbour = ancestor.HasLargerorEqual(targetmatid);
        if(!neighbour) DebugStop();
        if(!neighbour.Element()->Reference())
        {
            neighbour = ancestor.HasLarger(targetmatid);
        }
        while(!neighbour.Element()->Reference())
        {
            ancestor.SetCurrent(neighbour);
            neighbour = ancestor.HasLarger(targetmatid);
            if(!neighbour) break;
        }
        if(!neighbour) DebugStop();
        TPZCompEl *cellarge = neighbour.Element()->Reference();
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        TPZInterpolatedElement *intellarge = dynamic_cast<TPZInterpolatedElement *>(cellarge);
        intel->RestrainSide(gelside.Side(), intellarge, neighbour.Side());
    }
}

/// Set the Lagrange levels of the connects of the different meshes
void TPZMHMixedMeshControl::SetLagrangeMultiplierLevels()
{
    // MHMHDiv - there is "always" a lagrange multiplier mesh at the element level and
    // MHM level

    //    internal flux 0
    //    1, 2 or 3 displacement connects 1
    //    distributed element flux 2
    //    delayed internal pressure 3
    //    rotation 4
    //    boundary flux 5
    //    average element pressure 6
    //    delayed average element pressure except 1 7
    //    distributed MHM flux 8
    //    last delayed average element pressure 9
    //    MHM skeleton flux 10
    //    average MHM pressure 11
    //    MHM skeleton pressure 12

    int dim = fGMesh->Dimension();
    {
        int npress_delay = 0;
        switch(fProblemType)
        {
            case EScalar:
                npress_delay = 1;
                break;
            case EElasticity2D:
                npress_delay = 2;
                break;
            case EElasticity3D:
                npress_delay = 3;
                break;
            default:
                DebugStop();
        }
        int64_t nel = fPressureFineMesh->NElements();
        for (int64_t el = 0; el<nel; el++) {
            TPZCompEl *cel = fPressureFineMesh->Element(el);
            int subdomain = WhichSubdomain(cel);
            TPZGeoEl *gel = cel->Reference();
            int matid = gel->MaterialId();
            if(fMaterialIds.find(matid) != fMaterialIds.end())
            {
                int count_delay = 0;
                int ncon = cel->NConnects();
                for (int ic=0; ic<ncon; ic++) {
                    TPZConnect &c = cel->Connect(ic);
                    if(count_delay < npress_delay){
                        c.SetLagrangeMultiplier(1);
                        count_delay++;
                    } else
                    {
                        c.SetLagrangeMultiplier(3);
                    }
                }
            } else if(matid == fPressureSkeletonMatId)
            {
                int ncon = cel->NConnects();
                for (int ic=0; ic<ncon; ic++) {
                    TPZConnect &c = cel->Connect(ic);
                    c.SetLagrangeMultiplier(12);
                }
            } else {
                DebugStop();
            }
        }
    }
    // set the lagrange levels of the flux mesh
    {
        int64_t nel = fFluxMesh->NElements();
        for (int64_t el = 0; el<nel; el++) {
            TPZCompEl *cel = fFluxMesh->Element(el);
            int subdomain = WhichSubdomain(cel);
            TPZGeoEl *gel = cel->Reference();
            int matid = gel->MaterialId();
            if(fMaterialBCIds.find(matid) != fMaterialBCIds.end() || matid == fSkeletonMatId)
            {
                int ncon = cel->NConnects();
                for (int ic=0; ic<ncon; ic++) {
                    TPZConnect &c = cel->Connect(ic);
                    c.SetLagrangeMultiplier(10);
                }
            }
            else if (fMaterialIds.find(matid) != fMaterialIds.end())
            {
                int ncon = cel->NConnects();
                for (int ic=0; ic<ncon-1; ic++) {
                    TPZConnect &c = cel->Connect(ic);
                    c.SetLagrangeMultiplier(5);
                }
                // internal flux
                TPZConnect &c = cel->Connect(ncon-1);
                c.SetLagrangeMultiplier(0);
            }
            else if(matid == fWrapMatId)
            {
                
            }
            else
            {
                DebugStop();
            }
        }
    }
    if(fNState > 1)
    {
        int64_t nel = fRotationMesh->NElements();
        for (int64_t el = 0; el<nel; el++) {
            TPZCompEl *cel = fRotationMesh->Element(el);
            int ncon = cel->NConnects();
            for (int ic=0; ic<ncon; ic++) {
                TPZConnect &c = cel->Connect(ic);
                c.SetLagrangeMultiplier(4);
            }
        }
    }
    if(fLagrangeAveragePressure)
    {
        int64_t nel = fCMeshConstantPressure->NElements();
        std::map<int,int> delayed_pressure;
        for (int64_t el = 0; el<nel; el++) {
            TPZCompEl *cel = fCMeshConstantPressure->Element(el);
            int ncon = cel->NConnects();
            for (int ic=0; ic<ncon; ic++) {
                TPZConnect &c = cel->Connect(ic);
                c.SetLagrangeMultiplier(11);
            }
            cel = fCMeshLagrange->Element(el);
            for (int ic=0; ic<ncon; ic++) {
                TPZConnect &c = cel->Connect(ic);
                c.SetLagrangeMultiplier(8);
            }
        }
    }
    else
    {
        DebugStop();
    }
    {
        int64_t nel = fAverageSolution->NElements();
        if(!nel) DebugStop();
        std::map<int,int> delayed_pressure;
        for (int64_t el = 0; el<nel; el++) {
            TPZCompEl *cel = fAverageSolution->Element(el);
            int subdomain = WhichSubdomain(cel);
            int ncon = cel->NConnects();
            for (int ic=0; ic<ncon; ic++) {
                TPZConnect &c = cel->Connect(ic);
                if(delayed_pressure[subdomain] == 0)
                {
                    c.SetLagrangeMultiplier(9);
                    delayed_pressure[subdomain]++;
                }
                else
                {
                    c.SetLagrangeMultiplier(6);
                }
            }
            cel = fDistributedFlux->Element(el);
            for (int ic=0; ic<ncon; ic++) {
                TPZConnect &c = cel->Connect(ic);
                c.SetLagrangeMultiplier(2);
            }
        }
    }

}
