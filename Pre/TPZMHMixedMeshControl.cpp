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
#include "TPZSimpleTimer.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.mhmixedmeshcontrol");
#endif

/*
TPZMHMixedMeshControl::TPZMHMixedMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh, std::set<int64_t> &coarseindices) : TPZMHMeshControl(gmesh,coarseindices)
{
}
*/

TPZMHMixedMeshControl::TPZMHMixedMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh, TPZVec<int64_t> &coarseindices) : TPZMHMeshControl(gmesh,coarseindices)
{
}

TPZMHMixedMeshControl::TPZMHMixedMeshControl(int dimension) : TPZMHMeshControl(dimension)
{
}

TPZMHMixedMeshControl::TPZMHMixedMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh) : TPZMHMeshControl(gmesh)
{
}

TPZMHMixedMeshControl::~TPZMHMixedMeshControl()
{
}

/// Insert Boundary condition objects that do not perform any actual computation
void TPZMHMixedMeshControl::InsertPeriferalMaterialObjects()
{
    int matid = *fMaterialIds.begin();
    TPZMaterial *mat = fCMesh->FindMaterial(matid);
    if (!mat) {
        DebugStop();
    }
    if(fCMesh->MaterialVec().find(fSkeletonMatId) != fCMesh->MaterialVec().end())
    {
        std::cout << "Peripheral material inserted twice " << __PRETTY_FUNCTION__ << std::endl;
        DebugStop();
    }
    int nstate = mat->NStateVariables();
    int dim = fGMesh->Dimension();
    
    auto *nullmat = new TPZNullMaterialCS<STATE>(fSkeletonMatId);
    nullmat->SetDimension(mat->Dimension()-1);
    nullmat->SetNStateVariables(nstate);
    fCMesh->InsertMaterialObject(nullmat);

    // insert in flux mesh
    {
        for(auto matid : fMaterialIds)
        {
            TPZNullMaterial<> *mat = new TPZNullMaterial<>(matid,dim,nstate);
            fFluxMesh->InsertMaterialObject(mat);
        }
        for(auto matid : fMaterialBCIds)
        {
            TPZNullMaterial<> *mat = new TPZNullMaterial<>(matid,dim-1,nstate);
            fFluxMesh->InsertMaterialObject(mat);
        }
        TPZNullMaterial<> *mat = new TPZNullMaterial<>(fSkeletonMatId,dim-1,nstate);
        fFluxMesh->InsertMaterialObject(mat);
    }
    // insert in pressure mesh
    {
        for(auto matid : fMaterialIds)
        {
            TPZNullMaterial<> *mat = new TPZNullMaterial<>(matid,dim,nstate);
            fPressureFineMesh->InsertMaterialObject(mat);
        }
    }
    // insert in distributed flux (per element) mesh
    // insert in average solution (per element) mesh
    if(fHDivFamily == HDivFamily::EHDivStandard)
    {
    // insert in distributed flux per domain mesh
    // insert in average solution per domain mesh
        if(!fElementFluxMesh) fElementFluxMesh = new TPZCompMesh(fGMesh);
        if(!fElementAverageMesh) fElementAverageMesh = new TPZCompMesh(fGMesh);
        int nstate = 1;
        if(fProblemType == EElasticity2D) nstate = 3;
        if(fProblemType == EElasticity3D) nstate = 6;
        for(auto matid : fMaterialIds)
        {
            TPZNullMaterial<> *mat = new TPZNullMaterial<>(matid,dim,nstate);
            fElementFluxMesh->InsertMaterialObject(mat);
            TPZNullMaterial<> *mat2 = new TPZNullMaterial<>(matid,dim,nstate);
            fElementAverageMesh->InsertMaterialObject(mat2);
        }
    }
    
    {
        if(!fCMeshDomainFlux) fCMeshDomainFlux = new TPZCompMesh(fGMesh);
        if(!fCMeshDomainPressure) fCMeshDomainPressure = new TPZCompMesh(fGMesh);
        // insert material objects in the domain flux and average meshes
        for(auto matid : fMaterialIds)
        {
            int nstate = 1;
            if(fProblemType == EElasticity2D) nstate = 3;
            if(fProblemType == EElasticity3D) nstate = 6;
            TPZNullMaterial<> *mat = new TPZNullMaterial<>(matid,dim,nstate);
            fCMeshDomainFlux->InsertMaterialObject(mat);
            TPZNullMaterial<> *mat2 = new TPZNullMaterial<>(matid,dim,nstate);
            fCMeshDomainPressure->InsertMaterialObject(mat2);
        }
    }
    
    // insert material in the rotational mesh
    if(fProblemType == EElasticity2D || fProblemType == EElasticity3D)
    {
        if(!fRotationMesh) fRotationMesh = new TPZCompMesh(fGMesh);
        int nstate = 1;
        if(fProblemType == EElasticity3D) nstate = 3;
        for(auto matid : fMaterialIds)
        {
            TPZNullMaterial<> *mat = new TPZNullMaterial<>(matid,dim,nstate);
            fRotationMesh->InsertMaterialObject(mat);
        }
    }
}



/// Create all data structures for the computational mesh
void TPZMHMixedMeshControl::BuildComputationalMesh(bool usersubstructure)
{
    if (fpOrderInternal == 0 || fpOrderSkeleton == 0) {
        DebugStop();
    }
    
    std::cout << "\n----------- Creating atomic and multiphysics meshes -----------" << std::endl;
    TPZSimpleTimer cmeshestimes;

    
    // insert the peripheral material objects and create the auxiliary computational meshes if necessary
    InsertPeriferalMaterialObjects();
    CreateFluxMesh();
//    ChangeInternalOrder(fFluxMesh.operator->(), 2);
    
    
#ifdef PZDEBUG
    if (fFluxMesh->Dimension() != fGMesh->Dimension()) {
        DebugStop();
    }
#endif
    
    CreatePressureMHMMesh();
    if(fProblemType == EElasticity2D || fProblemType == EElasticity3D)
    {
        CreateRotationMesh();
    }
    
    if(fHDivFamily == HDivFamily::EHDivStandard || fProblemType == EElasticity2D || fProblemType == EElasticity3D)
    {
        CreateAverageElementMeshes();
    }
    
    CreateLagrangeMultiplierMesh();
    
    InitializeLagrangeLevels();
    
    CreateHDivPressureMHMMesh();
    
    std::cout << "Create cmeshes time: " << cmeshestimes.ReturnTimeDouble()/1000 << " seconds" << std::endl;
          
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
    
    std::cout << "\n----------- Before Hide the elements -----------" << std::endl;
    std::cout << "Number of equations in cmesh = " << fCMesh->NEquations() << std::endl;
    std::cout << "Number of rows in sol vec = " << fCMesh->Solution().Rows() << std::endl;
    std::cout << "\n----------- Starting Hiding the elements -----------" << std::endl;
    TPZSimpleTimer hidetime;
    if (usersubstructure) {
        HideTheElements();
    }
    std::cout << "HideTeElements time: " << hidetime.ReturnTimeDouble()/1000 << " seconds" << std::endl;
    std::cout << "Number of equations in cmesh = " << fCMesh->NEquations() << std::endl;
    std::cout << "Number of rows in sol vec = " << fCMesh->Solution().Rows() << std::endl;
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


void TPZMHMixedMeshControl::CreateFluxMesh()
{
    if(!fFluxMesh) DebugStop();
    TPZCompMesh * cmeshHDiv = fFluxMesh.operator->();
    cmeshHDiv->SetName("FluxMesh");

    CreateFluxElements();
    
#ifdef PZDEBUG
    if(0)
    {
        std::ofstream out("FluxmeshBeforeSkeleton.txt");
        fFluxMesh->Print(out);
    }
#endif
    
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



void TPZMHMixedMeshControl::CreatePressureMHMMesh()
{
    TPZGeoMesh * gmesh = fGMesh.operator->();
    gmesh->ResetReference();
    if(!fPressureFineMesh) DebugStop();
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
    int meshdim = gmesh->Dimension();
    int porder = fpOrderInternal;
    TPZCompMesh * cmeshRotation = fRotationMesh.operator->();
    gmesh->ResetReference();
    cmeshRotation->SetName("RotationMesh");
    cmeshRotation->SetDimModel(gmesh->Dimension());
    cmeshRotation->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
    cmeshRotation->SetDefaultOrder(porder);
    
    cmeshRotation->AutoBuild();
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
            DebugStop();
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




void TPZMHMixedMeshControl::CreateHDivPressureMHMMesh()
{
    TPZStack<TPZCompMesh *,3 > cmeshes;
    cmeshes.Push(fFluxMesh.operator->());
    cmeshes.Push(fPressureFineMesh.operator->());
    
#ifdef PZDEBUG
    if(1)
    {
        std::ofstream out("PressureMesh.txt");
        cmeshes[1] ->Print(out);
        
        std::ofstream out2("FluxMesh.txt");
        cmeshes[0] ->Print(out2);
    }
#endif
    
    if(fNState > 1) {
        cmeshes.Push(fRotationMesh.operator->());
#ifdef PZDEBUG
        if(1){
            std::ofstream out("RotationMesh.txt");
            fRotationMesh->Print(out);
        }
#endif
    }
    if(fHDivFamily == HDivFamily::EHDivStandard)
    {
#ifdef PZDEBUG
        if(1){
            std::ofstream out("ElFluxMesh.txt");
            fElementFluxMesh->Print(out);
            std::ofstream out2("ElPressureMesh.txt");
            fElementAverageMesh->Print(out2);
        }
#endif
        cmeshes.Push(fElementFluxMesh.operator->());
        cmeshes.Push(fElementAverageMesh.operator->());
    }
    
    // cmeshes.Push(fCMeshDomainFlux.operator->());
    // cmeshes.Push(fCMeshDomainPressure.operator->());
#ifdef PZDEBUG
    if(1){
        std::ofstream out("DomainFluxMesh.txt");
        fCMeshDomainFlux->Print(out);
        std::ofstream out2("DomainPressureMesh.txt");
        fCMeshDomainPressure->Print(out2);
    }
#endif

    TPZGeoMesh *gmesh = cmeshes[0]->Reference();
    if(!gmesh)
    {
        std::cout<< "Geometric mesh doesn't exist" << std::endl;
        DebugStop();
    }
    int dim = gmesh->Dimension();
    gmesh->ResetReference();
    // Multiphysics mesh
    TPZMultiphysicsCompMesh * MixedFluxPressureCmesh = fCMesh.operator->();
    MixedFluxPressureCmesh->SetDimModel(dim);
    MixedFluxPressureCmesh->SetAllCreateFunctionsMultiphysicElem();
    
    fCMesh->BuildMultiphysicsSpace(cmeshes);
    TPZManVector<TPZCompMesh * ,7> meshvector;
    
#ifdef PZDEBUG
    if(0)
    {
        std::ofstream out2("gmesh.txt");
        gmesh->Print(out2);
    }
#endif
    
    meshvector = cmeshes;

    // populate the connect to subdomain data structure for the multiphysics mesh
    JoinSubdomains(meshvector, MixedFluxPressureCmesh);
    MixedFluxPressureCmesh->ComputeNodElCon();

    // Transferindo para a multifisica
//    TPZBuildMultiphysicsMesh::AddElements(meshvector, MixedFluxPressureCmesh);
//    TPZBuildMultiphysicsMesh::AddConnects(meshvector, MixedFluxPressureCmesh);
    // copy the solution of the atomic meshes to the multiphysics mesh

#ifdef PZDEBUG
    if(0)
    {
        std::ofstream file("cmeshmphys.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(MixedFluxPressureCmesh, file,true);
        std::ofstream out("cmeshmphys.txt");
        MixedFluxPressureCmesh->Print(out);
    }
#endif
    
#ifdef PZDEBUG
    if(1)
    {
        std::ofstream out("MultiphysicsMesh.txt");
        MixedFluxPressureCmesh->Print(out);
    }
#endif
    
    return;
    
}

void TPZMHMixedMeshControl::HideTheElements()
{
//    if (fHybridize) {
//        KeepOneLagrangian = false;
//    }
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
    int64_t count = 0;
    for(auto &it : ElementGroups) {
        int64_t index;
        int KeepOneLagrangian = LagrangeLevels::DomainPressure;
        TPZCompMeshTools::PutinSubmeshes(fCMesh.operator->(), it.second, index, KeepOneLagrangian);
        submeshindices[it.first] = index;
    }
//    TPZCompMeshTools::PutinSubmeshes(fCMesh.operator->(), ElementGroups, submeshindices, KeepOneLagrangian);
    std::cout << "After putting in substructures\n";
    fMHMtoSubCMesh = submeshindices;
    fCMesh->ComputeNodElCon();
    // if there are constrained connects that have been transferred
    {
        int64_t ncon = fCMesh->NConnects();
        for(int64_t ic = 0; ic<ncon; ic++)
        {
            TPZConnect &c = fCMesh->ConnectVec()[ic];
            if(c.NElConnected() == 0 && c.HasDependency())
            {
                c.RemoveDepend();
            }
        }
    }
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
    if (fCMeshDomainFlux)
    {
        out << "******************* LAGRANGE MULTIPLIER MESH *****************\n";
        fCMeshDomainFlux->Print(out);
    }
    
    /// computational mesh to represent the constant states
    if (fCMeshDomainPressure)
    {
        out << "******************* CONSTANTE PRESSURE MESH *****************\n";
        fCMeshDomainPressure->Print(out);
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
    
}

/*
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

    DebugStop();

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
*/
    
// create the elements domain per domain with approximation spaces disconnected from each other
void TPZMHMixedMeshControl::CreateFluxElements() {
    TPZGeoMesh *gmesh = fGMesh.operator->();
    int meshdim = gmesh->Dimension();
    TPZCompMesh * cmeshHDiv = fFluxMesh.operator->();
    gmesh->ResetReference();
    cmeshHDiv->LoadReferences();
    cmeshHDiv->SetDimModel(meshdim);
    cmeshHDiv->ApproxSpace().SetAllCreateFunctionsHDiv(meshdim);
    cmeshHDiv->SetDefaultOrder(fpOrderInternal);
    
    std::set<int> matids;
    std::merge(fMaterialIds.begin(), fMaterialIds.end(),
                    fMaterialBCIds.begin(), fMaterialBCIds.end(),
                    std::inserter(matids, matids.begin()));
    matids.insert(fSkeletonMatId);
    cmeshHDiv->AutoBuild(matids);
    fConnectToSubDomainIdentifier[cmeshHDiv].Resize(cmeshHDiv->NConnects(), -1);
    //Criar elementos computacionais malha MHM
    TPZGeoEl *gel = NULL;
    TPZGeoEl *gsubel = NULL;
    int64_t nel = cmeshHDiv->NElements();
    for(int64_t el = 0; el<nel; el++)
    {
        TPZCompEl *cel = cmeshHDiv->Element(el);
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        TPZGeoEl *gel = cel->Reference();
        int eldim = gel->Dimension();
        if(eldim < meshdim) {
            int domain = fGeoToMHMDomain[gel->Index()];
            SetSubdomain(cel, domain);
            if(gel->MaterialId() == fSkeletonMatId)
            {
                intel->PRefine(fpOrderSkeleton);
                intel->SetPreferredOrder(fpOrderSkeleton);
            }
            continue;
        }
        int nfaces = gel->NSides(meshdim-1);
        int nsides = gel->NSides();
        for (int is = nsides-nfaces-1; is<nsides-1; is++) {
            int sidedomain = fGeoToMHMDomain[gel->Index()];
            if(sidedomain == -1) DebugStop();
            int64_t cindex = intel->SideConnectIndex(0, is);
            TPZGeoElSide gelside(gel,is);
            TPZGeoElSide neighbour = gelside.Neighbour();
            while(neighbour != gelside)
            {
                int neighdomain = fGeoToMHMDomain[neighbour.Element()->Index()];
                if(neighdomain != -1 && neighdomain != sidedomain)
                {
                    sidedomain = -1;
                }
                neighbour = neighbour.Neighbour();
            }
            if(sidedomain != -1)
            {
                SetSubdomain(cmeshHDiv, cindex, sidedomain);
            }
        }
        if(fHdivmaismais > 0)
        {
            intel->ForceSideOrder(nsides-1,fpOrderInternal+fHdivmaismais);
        }
    }
    
    fFluxMesh->ExpandSolution();
    
}

/// Create the meshes which represent element distributed flux and element average
void TPZMHMixedMeshControl::CreateAverageElementMeshes()
{
    if(fGMesh->Reference()) fGMesh->ResetReference();
    {
        fElementFluxMesh->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
        fElementFluxMesh->SetDefaultOrder(0);
        fElementFluxMesh->AutoBuild();
        int64_t nel = fElementFluxMesh->NElements();
        for (int64_t el = 0; el<nel; el++) {
            TPZCompEl *cel = fElementFluxMesh->Element(el);
            TPZGeoEl *gel = cel->Reference();
            int64_t index = gel->Index();
            int domain = fGeoToMHMDomain[index];
            SetSubdomain(cel, domain);
        }
    }
    fGMesh->ResetReference();
    {
        fElementAverageMesh->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
        fElementAverageMesh->SetDefaultOrder(0);
        fElementAverageMesh->AutoBuild();
        int64_t nel = fElementAverageMesh->NElements();
        for (int64_t el = 0; el<nel; el++) {
            TPZCompEl *cel = fElementAverageMesh->Element(el);
            TPZGeoEl *gel = cel->Reference();
            int64_t index = gel->Index();
            int domain = fGeoToMHMDomain[index];
            SetSubdomain(cel, domain);
        }
    }
    fGMesh->ResetReference();
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
        char keeplagrange = LagrangeLevels::Elpressure;
        TPZCompMeshTools::CondenseElements(subcmesh, keeplagrange, true);
//        TPZCompMeshTools::CondenseElements(subcmesh, keeplagrange);
        subcmesh->CleanUpUnconnectedNodes();
        int numthreads = 16;
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
//        subcmesh->SetAnalysisFStruct(numthreads);
//        subcmesh->SetAnalysisSkyline(numthreads, preconditioned, guiInterface);
        subcmesh->SetAnalysisSparse(numthreads);
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

// set the Lagrange levels of the connects for the atomic meshes
void TPZMHMixedMeshControl::InitializeLagrangeLevels()
{
    auto meshvec = GetMeshes();
    auto fluxmesh = meshvec[0];
    int dim = fluxmesh->Dimension();
    // if(fProblemType != EElasticity3D || meshvec.size() != 7)
    // {
    //     std::cout << "Configuration not tested, calling DebugStop for you to verify the connect levels\n";
    //     DebugStop();
    // }
    {
        int64_t nel = fluxmesh->NElements();
        for(int64_t el = 0; el<nel; el++)
        {
            TPZCompEl *cel = fluxmesh->Element(el);
            if(!cel) continue;
            TPZGeoEl *gel = cel->Reference();
            int eldim = gel->Dimension();
            if(eldim == dim)
            {
                int nc = cel->NConnects();
                for(int ic = 0; ic<nc-1; ic++) cel->Connect(ic).SetLagrangeMultiplier(5);
                cel->Connect(nc-1).SetLagrangeMultiplier(0);
            }
            else if(gel->MaterialId() == fSkeletonMatId)
            {
                cel->Connect(0).SetLagrangeMultiplier(9);
            }
        }
    }
    {
        auto pressuremesh = meshvec[1];
        int64_t nel = pressuremesh->NElements();
        for(int64_t el = 0; el<nel; el++)
        {
            TPZCompEl *cel = pressuremesh->Element(el);
            if(!cel || cel->Reference()->Dimension() != dim) continue;
            int nc = cel->NConnects();
            int nfirstconnects = 3;
            for(int ic = 0; ic<nfirstconnects; ic++) cel->Connect(ic).SetLagrangeMultiplier(1);
            for(int ic = nfirstconnects; ic<nc; ic++) cel->Connect(ic).SetLagrangeMultiplier(3);
//            cel->Connect(nc-1).SetLagrangeMultiplier(5);
        }
    }
    if(fRotationMesh)
    {
        int64_t nel = fRotationMesh->NElements();
        for(int64_t el = 0; el<nel; el++)
        {
            TPZCompEl *cel = fRotationMesh->Element(el);
            if(!cel || cel->Reference()->Dimension() != dim) continue;
            int nc = cel->NConnects();
            if(nc != 1) DebugStop();
            cel->Connect(0).SetLagrangeMultiplier(4);
        }
    }
    {
        int64_t nel = fElementFluxMesh->NElements();
        for(int64_t el = 0; el<nel; el++)
        {
            TPZCompEl *cel = fElementFluxMesh->Element(el);
            if(!cel || cel->Reference()->Dimension() != dim) continue;
            int nc = cel->NConnects();
            if(nc != 1) DebugStop();
            cel->Connect(0).SetLagrangeMultiplier(2);
        }
    }
    {
        int64_t nel = fElementAverageMesh->NElements();
        std::set<int64_t> subdomains;
        for(int64_t el = 0; el<nel; el++)
        {
            TPZCompEl *cel = fElementAverageMesh->Element(el);
            if(!cel || cel->Reference()->Dimension() != dim) continue;
            int nc = cel->NConnects();
            TPZGeoEl *gel = cel->Reference();
            int64_t domain = fGeoToMHMDomain[gel->Index()];
            if(nc != 1) DebugStop();
            cel->Connect(0).SetLagrangeMultiplier(6);
            // this should be one element per subdomain
            // one element por subdomain has level 8
            if(subdomains.find(domain) ==  subdomains.end())
            {
                subdomains.insert(domain);
                cel->Connect(0).SetLagrangeMultiplier(8);
            }
        }
    }
    {
        int64_t nel = fCMeshDomainFlux->NElements();
        for(int64_t el = 0; el<nel; el++)
        {
            TPZCompEl *cel = fCMeshDomainFlux->Element(el);
            if(!cel || cel->Reference()->Dimension() != dim) continue;
            int nc = cel->NConnects();
            if(nc != 1) DebugStop();
            cel->Connect(0).SetLagrangeMultiplier(7);
        }
    }
    {
        int64_t nel = fCMeshDomainPressure->NElements();
        for(int64_t el = 0; el<nel; el++)
        {
            TPZCompEl *cel = fCMeshDomainPressure->Element(el);
            if(!cel || cel->Reference()->Dimension() != dim) continue;
            int nc = cel->NConnects();
            if(nc != 1) DebugStop();
            cel->Connect(0).SetLagrangeMultiplier(10);
        }
    }

}

void TPZMHMixedMeshControl::ChangeInternalOrder(TPZCompMesh *cmesh, int pOrder) {
    int64_t nel = cmesh->NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        if (!cel) continue;

        TPZGeoEl *gel = cel->Reference();
        if (gel->Dimension() != cmesh->Dimension()) {
            continue;
        }
        int nc = cel->NConnects();
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cel);
        if (!intel) DebugStop();

        intel->ForceSideOrder(gel->NSides() - 1, pOrder);
    }
    cmesh->ExpandSolution();
}
