//
//  TPZMHMixedMeshControl.cpp
//  PZ
//
//  Created by Philippe Devloo on 09/10/16.
//
//

#include "TPZMHMixedMeshControl.h"
#include "TPZBndCond.h"
#include "TPZNullMaterial.h"
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
#include "TPZNullMaterialCS.h"
#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.mhmixedmeshcontrol");
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
    fFluxMesh = new TPZCompMesh(gmesh);
    fFluxMesh->SetDimModel(gmesh->Dimension());
    fCMesh = new TPZMultiphysicsCompMesh(gmesh);
    fCMesh->SetDimModel(gmesh->Dimension());
}

TPZMHMixedMeshControl::TPZMHMixedMeshControl(int dimension) : TPZMHMeshControl(dimension)
{
    fFluxMesh = new TPZCompMesh(fGMesh);
    fFluxMesh->SetDimModel(fGMesh->Dimension());
    fCMesh = new TPZMultiphysicsCompMesh(fGMesh);
    fCMesh->SetDimModel(fGMesh->Dimension());
}

TPZMHMixedMeshControl::TPZMHMixedMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh) : TPZMHMeshControl(gmesh)
{
    fFluxMesh = new TPZCompMesh(gmesh);
    fFluxMesh->SetDimModel(gmesh->Dimension());
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
    int matid = *fMaterialIds.begin();
    TPZMaterial *mat = fCMesh->FindMaterial(matid);
    if (!mat) {
        DebugStop();
    }
    int nstate = mat->NStateVariables();
    TPZFNMatrix<1,STATE> val1(nstate,nstate,0.), val2Flux(nstate,1,0.);
    int typePressure = 0;
    
    if(fCMesh->MaterialVec().find(fSkeletonMatId) != fCMesh->MaterialVec().end())
    {
        std::cout << "Peripheral material inserted twice " << __PRETTY_FUNCTION__ << std::endl;
        return;
    }
    
    auto *nullmat = new TPZNullMaterialCS<STATE>(fSkeletonMatId);
    nullmat->SetDimension(mat->Dimension()-1);
    nullmat->SetNStateVariables(nstate);
    fCMesh->InsertMaterialObject(nullmat);
//    TPZBndCond * bcFlux = mat->CreateBC(mat, fSkeletonMatId, typePressure, val1, val2Flux);
    //    bcN->SetForcingFunction(0,force);
//    fCMesh->InsertMaterialObject(bcFlux);
    
    


}



/// Create all data structures for the computational mesh
void TPZMHMixedMeshControl::BuildComputationalMesh(bool usersubstructure)
{
    if (fpOrderInternal == 0 || fpOrderSkeleton == 0) {
        DebugStop();
    }
    InsertPeriferalMaterialObjects();
    CreateHDivMHMMesh();
    
    InsertPeriferalPressureMaterialObjects();
    if(fNState > 1){
        fRotationMesh = new TPZCompMesh(fGMesh);
        InsertPeriferalRotationMaterialObjects();
    }
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
    
    CreateHDivPressureMHMMesh();
    std::cout << "Total number of equations " << fCMesh->Solution().Rows() << std::endl;
    fGlobalSystemWithLocalCondensationSize = fCMesh->NEquations();
    fGlobalSystemSize = fCMesh->Solution().Rows();

    fCMesh->ComputeNodElCon();
#ifdef PZDEBUG2
    {
        std::ofstream out("Friendly.txt");
        PrintFriendly(out);
    }
    CheckMeshConsistency();
#endif
    
    if (usersubstructure) {
        HideTheElements();
    }
    fNumeq = fCMesh->NEquations();
    
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
    CreateSkeleton();
    
    if (fHdivmaismais) {
        int64_t nel = cmeshHDiv->ElementVec().NElements();
        for (int64_t el = 0; el < nel; el++) {
            TPZCompEl *cel = cmeshHDiv->ElementVec()[el];
            TPZGeoEl *gel = cel->Reference();
            if (gel->Dimension() == cmeshHDiv->Dimension()) {
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
    if(0)
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
    gmesh->ResetReference();
    cmeshHDiv->LoadReferences();
    cmeshHDiv->SetDimModel(meshdim);
    cmeshHDiv->ApproxSpace().SetAllCreateFunctionsHDiv(meshdim);
    cmeshHDiv->SetDefaultOrder(fpOrderInternal);
    if(fMaterialIds.size() == 0) DebugStop();
    TPZMaterialT<STATE> *matl2;
    for (auto item:fMaterialIds) {
        auto *mat = new TPZNullMaterial<STATE>(item);
        matl2 = mat;
        mat->SetNStateVariables(fNState);
        mat->SetDimension(meshdim);
        cmeshHDiv->InsertMaterialObject(mat);
    }
    for (auto item:fMaterialBCIds) {
        TPZFNMatrix<1,STATE> val1(fNState,fNState,0.);
        TPZManVector<STATE,1>val2(fNState,0.);
        TPZBndCondT<STATE> *bc = matl2->CreateBC(matl2, item, 0, val1, val2);
        cmeshHDiv->InsertMaterialObject(bc);
    }
    {
        auto *nullmat = new TPZNullMaterial<STATE>(fSkeletonMatId);
        nullmat->SetDimension(meshdim-1);
        nullmat->SetNStateVariables(fNState);
        cmeshHDiv->InsertMaterialObject(nullmat);
        if(fSecondSkeletonMatId != 0)
        {
            auto *nullmat = new TPZNullMaterial<STATE>(fSecondSkeletonMatId);
            nullmat->SetDimension(meshdim-1);
            nullmat->SetNStateVariables(fNState);
            cmeshHDiv->InsertMaterialObject(nullmat);
        }
    }

}



void TPZMHMixedMeshControl::CreatePressureMHMMesh()
{
    TPZGeoMesh * gmesh = fGMesh.operator->();
    gmesh->ResetReference();

    // the pressure mesh should be empty when calling this method
    int64_t nskeletonconnects = fPressureFineMesh->NConnects();
    if(nskeletonconnects != 0){
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

    
#ifdef PZDEBUG
    // a very strange check!! Why does material id 1 need to be volumetric?
    {
        int64_t nel = fGMesh->NElements();
        for (int64_t el = 0; el<nel; el++) {
            TPZGeoEl *gel = fGMesh->Element(el);
            if (gel && gel->MaterialId() == 1) {
                if (gel->Dimension() != fGMesh->Dimension()) {
                    DebugStop();
                }
            }
        }
    }
#endif
    
    // the lagrange multiplier level is set to one
    int64_t nc = cmeshPressure->NConnects();
    for (int64_t ic=nskeletonconnects; ic<nc; ic++) {
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
    
    if(0){
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
            auto *matl2 = new TPZNullMaterial<STATE>((matid));
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
        auto *mathybrid = new TPZNullMaterial<STATE>(fPressureSkeletonMatId);
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
            auto *matl2 = new TPZNullMaterial<STATE>((matid));
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



void TPZMHMixedMeshControl::CreateHDivPressureMHMMesh()
{
    TPZManVector<TPZCompMesh *,3 > cmeshes(2);
    cmeshes[0] = fFluxMesh.operator->();
    cmeshes[1] = fPressureFineMesh.operator->();
    
    if(0)
    {
        std::ofstream out("PressureMesh_MultiPhis.txt");
        cmeshes[1] ->Print(out);
        
        std::ofstream out2("FluxMesh_MultiPhis.txt");
        cmeshes[0] ->Print(out2);
    }
    
    
    if(fNState > 1) {
        cmeshes.Resize(3);
        cmeshes[2] = fRotationMesh.operator->();
    }
    TPZGeoMesh *gmesh = cmeshes[0]->Reference();
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
    
    BuildMultiPhysicsMesh();
    TPZManVector<TPZCompMesh * ,3> meshvector;
    
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
    
    meshvector = cmeshes;

    // populate the connect to subdomain data structure for the multiphysics mesh
    JoinSubdomains(meshvector, MixedFluxPressureCmesh);
    
    // Transferindo para a multifisica
//    TPZBuildMultiphysicsMesh::AddElements(meshvector, MixedFluxPressureCmesh);
//    TPZBuildMultiphysicsMesh::AddConnects(meshvector, MixedFluxPressureCmesh);
    // copy the solution of the atomic meshes to the multiphysics mesh
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector, MixedFluxPressureCmesh);

    std::pair<int,int> skelmatid(fSkeletonMatId,fSecondSkeletonMatId);
    CreateMultiPhysicsInterfaceElements(fGMesh->Dimension()-1);
    CreateMultiPhysicsInterfaceElements(fGMesh->Dimension()-2);
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
    bool KeepOneLagrangian = true;
    if (fHybridize) {
        KeepOneLagrangian = false;
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
    TPZCompMeshTools::PutinSubmeshes(fCMesh.operator->(), ElementGroups, submeshindices, KeepOneLagrangian);
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
        TPZCompEl* celflux = fFluxMesh->CreateCompEl(skeleton2.CreatedElement());
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(celflux);
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
//        fFluxMesh->CreateCompEl(pressuregel.CreatedElement(), indexpressure);
        TPZCompEl* celpressure = fPressureFineMesh->CreateCompEl(pressuregel.CreatedElement());
        TPZInterpolatedElement *presel = dynamic_cast<TPZInterpolatedElement *>(celpressure);
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
    int64_t nel = fGMesh->NElements();

    if (fGeoToMHMDomain.size() != nel) DebugStop();

    TPZManVector<std::pair<int64_t, int64_t>> MHMOfEachGeoEl(nel);
    for (int i = 0; i < nel; i++) {
        MHMOfEachGeoEl[i] = std::make_pair(fGeoToMHMDomain[i], i);
    }

    std::sort(&MHMOfEachGeoEl[0], &MHMOfEachGeoEl[nel - 1] + 1);

    int64_t previousMHMDomain = -1;
    int64_t firstElemInMHMDomain = -1;
    for (int i = 0; i < MHMOfEachGeoEl.size(); i++) {
        int64_t MHMDomain = MHMOfEachGeoEl[i].first;
        int64_t elIndex = MHMOfEachGeoEl[i].second;

        if (MHMDomain == -1) continue;

        if (MHMDomain != previousMHMDomain) {
            if (previousMHMDomain != -1) {
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
        TPZCompEl* cel = cmeshHDiv->CreateCompEl(gel);
//        TPZCompEl *cel = cmeshHDiv->Element(index);
        // Associate the connects with the subdomain for the flux mesh
        SetSubdomain(cel, MHMDomain);
    }
    // Resets references of last MHM domain
    for (int j = firstElemInMHMDomain; j < MHMOfEachGeoEl.size(); j++) {
        fGMesh->Element(MHMOfEachGeoEl[j].second)->ResetReference();
    }

    fFluxMesh->ExpandSolution();
}

/// will create the elements on the skeleton
void TPZMHMixedMeshControl::CreateSkeleton()
{
    // comment this line or not to switch the type of skeleton elements
    int meshdim = fFluxMesh->Dimension();
    fFluxMesh->ApproxSpace().SetAllCreateFunctionsHDiv(meshdim);
    fGMesh->ResetReference();
    int order = fpOrderSkeleton;
    if (order <= 0) {
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
        // create an element to model the flux between subdomains
        TPZCompEl *cel = fFluxMesh->CreateCompEl(gel);
//        TPZCompEl *cel = fFluxMesh->ElementVec()[index];
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

#ifdef PZ_LOG
        if(logger.isDebugEnabled())
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


/*
TPZGeoMesh *gmesh = fGMesh.operator->();
int meshdim = gmesh->Dimension();
TPZCompMesh * cmeshHDiv = fFluxMesh.operator->();
gmesh->ResetReference();
cmeshHDiv->LoadReferences();
cmeshHDiv->SetDimModel(meshdim);
cmeshHDiv->ApproxSpace().SetAllCreateFunctionsHDiv(meshdim);
cmeshHDiv->SetDefaultOrder(fpOrderInternal);
TPZVecL2 *matl2 = new TPZVecL2(1);
cmeshHDiv->InsertMaterialObject(matl2);
matl2 = new TPZVecL2(2);
cmeshHDiv->InsertMaterialObject(matl2);
TPZFNMatrix<1,STATE> val1(1,1,0.),val2(1,1,0.);
TPZBndCond *bc = matl2->CreateBC(matl2, -1, 0, val1, val2);
cmeshHDiv->InsertMaterialObject(bc);
bc = matl2->CreateBC(matl2, -2, 0, val1, val2);
cmeshHDiv->InsertMaterialObject(bc);
bc = matl2->CreateBC(matl2, fSkeletonMatId, 0, val1, val2);
cmeshHDiv->InsertMaterialObject(bc);

int LagrangeMatIdLeft = 50;
int LagrangeMatIdRight = 51;
int nstate = 1;
TPZLagrangeMultiplier *matleft = new TPZLagrangeMultiplier(LagrangeMatIdLeft,meshdim,nstate);
TPZLagrangeMultiplier *matright = new TPZLagrangeMultiplier(LagrangeMatIdRight,meshdim,nstate);
matleft->SetMultiplier(-1.);
matright->SetMultiplier(-1.);
cmeshHDiv->InsertMaterialObject(matleft);
cmeshHDiv->InsertMaterialObject(matright);


std::set<int> materialids;
materialids.insert(1);
materialids.insert(-1);
materialids.insert(-2);
cmeshHDiv->AutoBuild(materialids);

cmeshHDiv->SetDefaultOrder(fpOrderSkeleton);
materialids.clear();
materialids.insert(fSkeletonMatId);
cmeshHDiv->AutoBuild(materialids);


// set the subdomain index of the connects for the flux elements
std::map<int64_t,int64_t>::iterator it;
for(it = fCoarseIndices.begin(); it != fCoarseIndices.end(); it++)
{
    TPZGeoEl *gel = fGMesh->Element(it->first);
    TPZStack<TPZCompElSide> connected;
    TPZGeoElSide gelside(gel,gel->NSides()-1);
    gelside.EqualorHigherCompElementList2(connected, 0, 0);
    int nst = connected.size();
    for (int64_t ist = 0; ist<nst; ist++) {
        TPZCompEl *cel = connected[ist].Element();
        SetSubdomain(cel, it->first);
    }
}
*/

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
            if (neighbour.Element()->MaterialId() == fSkeletonWrapMatId) {
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
//        celside.Element()->Print();
//        celstack[0].Element()->Print();
//        celstack[1].Element()->Print();
        TPZCompEl* cel1 = new TPZMultiphysicsInterfaceElement(*MixedFluxPressureCmesh,gbcleft.CreatedElement(),celstack[0],celside);
        TPZCompEl* cel2 = new TPZMultiphysicsInterfaceElement(*MixedFluxPressureCmesh,gbcright.CreatedElement(),celstack[1],celside);
#ifdef PZ_LOG
        if(logger.isDebugEnabled())
        {
            std::stringstream sout;
            sout << "Created an interface skeleton element index " << cel1->Index() << " between mphys " << celstack[0].Element()->Index() << " and pressure mphyx " << celside.Element()->Index()
            << std::endl;
            sout << "Created an interface skeleton element index " << cel2->Index() << " between mphys " << celstack[1].Element()->Index() << " and pressure mphyx " << celside.Element()->Index();
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
        
#ifdef PZ_LOG2
        if(logger.isDebugEnabled())
        {
            std::stringstream sout;
            subcmesh->Print(sout);
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        // TODO: increment nelconnected of exterior connects
        bool keeplagrange = true;
        TPZCompMeshTools::CreatedCondensedElements(subcmesh, keeplagrange);
        subcmesh->CleanUpUnconnectedNodes();
        int numthreads = 0;
        int preconditioned = 0;
#ifdef PZ_LOG2
        if(logger.isDebugEnabled())
        {
            std::stringstream sout;
            subcmesh->Print(sout);
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        TPZAutoPointer<TPZGuiInterface> guiInterface;
        subcmesh->SetAnalysisSkyline(numthreads, preconditioned, guiInterface);
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
#ifdef PZ_LOG
    if(logger.isDebugEnabled())
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

void TPZMHMixedMeshControl::AdjustBoundaryElements() {

    // store the entries of the map that correspond to boundary elements. these entries will be removed from the map
    std::set<int64_t> entries_to_erase;

    for (auto it = fInterfaces.cbegin(); it != fInterfaces.cend();) {
        int64_t elindex = it->first;
        TPZGeoEl *gel = fGMesh->Element(elindex);

        // filter boundary elements
        if (elindex == it->second.second) {

            // create boundary compels
            TPZCompEl* cel = fFluxMesh->CreateCompEl(gel);
            SetSubdomain(cel, -1);

            // set the side orientation of the boundary elements
            int side = gel->NSides() - 1;
            TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
            intel->SetSideOrient(side, 1);

            TPZGeoElSide gelside(gel);
            for (TPZGeoElSide neighbour = gelside.Neighbour(); neighbour != gelside; neighbour++) {
                TPZGeoEl *neigh_element = neighbour.Element();
                if (neigh_element->Dimension() != fGMesh->Dimension()) continue;
                SetSubdomain(cel, fGeoToMHMDomain[neigh_element->Index()]);
            }
            entries_to_erase.insert(it->first);
        }
        it++;
    }

    auto it = fInterfaces.cbegin();
    while (it != fInterfaces.cend()) {
        if (entries_to_erase.find(it->first) != entries_to_erase.cend()) {
            it = fInterfaces.erase(it);
        } else {
            ++it;
        }
    }
}

