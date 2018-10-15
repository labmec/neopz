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
#include "TPZMatLaplacian.h"
#include "TPZLagrangeMultiplier.h"


#include <iostream>
#include <sstream>
#include <iterator>
#include <numeric>

#include "pzsubcmesh.h"

#include "pzbuildmultiphysicsmesh.h"
#include "TPZCompMeshTools.h"
#include "pzinterpolationspace.h"
#include "pzintel.h"
#include "TPZInterfaceEl.h"
#include "TPZMultiphysicsInterfaceEl.h"

#include "TPZVTKGeoMesh.h"

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
    fFluxMesh = new TPZCompMesh(gmesh);
    fFluxMesh->SetDimModel(gmesh->Dimension());
    fCMesh = new TPZCompMesh(gmesh);
    fCMesh->SetDimModel(gmesh->Dimension());
}

TPZMHMixedMeshControl::TPZMHMixedMeshControl(int dimension) : TPZMHMeshControl(dimension)
{
    fFluxMesh = new TPZCompMesh(fGMesh);
    fFluxMesh->SetDimModel(fGMesh->Dimension());
    fCMesh = new TPZCompMesh(fGMesh);
    fCMesh->SetDimModel(fGMesh->Dimension());
}

TPZMHMixedMeshControl::TPZMHMixedMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh) : TPZMHMeshControl(gmesh)
{
    fFluxMesh = new TPZCompMesh(gmesh);
    fFluxMesh->SetDimModel(gmesh->Dimension());
    fCMesh = new TPZCompMesh(gmesh);
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
    TPZFNMatrix<1,STATE> val1(1,1,0.), val2Flux(1,1,0.);
    int typePressure = 0;
    
    TPZBndCond * bcFlux = mat->CreateBC(mat, fSkeletonMatId, typePressure, val1, val2Flux);
    //    bcN->SetForcingFunction(0,force);
    fCMesh->InsertMaterialObject(bcFlux);
    
    


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

#ifdef PZDEBUG
    if (fFluxMesh->Dimension() != fGMesh->Dimension()) {
        DebugStop();
    }
#endif
    
    CreatePressureMHMMesh();
    
    
    CreateHDivPressureMHMMesh();
    std::cout << "Total number of equations " << fCMesh->Solution().Rows() << std::endl;
    fGlobalSystemWithLocalCondensationSize = fCMesh->NEquations();
    fGlobalSystemSize = fCMesh->Solution().Rows();

    fCMesh->ComputeNodElCon();
#ifdef PZDEBUG
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

}


void TPZMHMixedMeshControl::CreateHDivMHMMesh()
{
    TPZCompMesh * cmeshHDiv = fFluxMesh.operator->();
    InsertPeriferalHdivMaterialObjects();
    CreateInternalFluxElements();
    
#ifdef PZDEBUG
    {
        fFluxMesh->ExpandSolution();
        std::ofstream out("FluxmeshBeforeSkeleton.txt");
        fFluxMesh->Print(out);
    }
#endif
    CreateSkeleton();


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
    fFluxMesh->ExpandSolution();
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
    TPZMaterial *matl2;
    for (auto item:fMaterialIds) {
        TPZVecL2 *mat = new TPZVecL2(item);
        matl2 = mat;
        cmeshHDiv->InsertMaterialObject(matl2);
    }
    for (auto item:fMaterialBCIds) {
        TPZFNMatrix<1,STATE> val1(1,1,0.),val2(1,1,0.);
        TPZBndCond *bc = matl2->CreateBC(matl2, item, 0, val1, val2);
        cmeshHDiv->InsertMaterialObject(bc);
    }
    {
        TPZFNMatrix<1,STATE> val1(1,1,0.),val2(1,1,0.);
        TPZBndCond *bc = matl2->CreateBC(matl2, fSkeletonMatId, 0, val1, val2);
        cmeshHDiv->InsertMaterialObject(bc);
        if(fSecondSkeletonMatId != 0)
        {
            bc = matl2->CreateBC(matl2, fSecondSkeletonMatId, 0, val1, val2);
            cmeshHDiv->InsertMaterialObject(bc);
        }
    }

}



void TPZMHMixedMeshControl::CreatePressureMHMMesh()
{
    TPZGeoMesh * gmesh = fGMesh.operator->();
    gmesh->ResetReference();
    
    int64_t nskeletonconnects = fPressureFineMesh->NConnects();
    int porder = fpOrderInternal;
    TPZCompMesh * cmeshPressure = fPressureFineMesh.operator->();
    gmesh->ResetReference();
    cmeshPressure->SetName("PressureMesh");
    cmeshPressure->SetDimModel(gmesh->Dimension());
    cmeshPressure->ApproxSpace().SetAllCreateFunctionsContinuous();
    cmeshPressure->ApproxSpace().CreateDisconnectedElements(true);
    cmeshPressure->SetDefaultOrder(porder);
    int meshdim = cmeshPressure->Dimension();
    std::set<int> matids;
    for (auto it:fMaterialIds) {
        TPZMaterial *mat = cmeshPressure->FindMaterial(it);
        if (mat && mat->Dimension() == meshdim) {
            matids.insert(it);
        }
    }
    cmeshPressure->AutoBuild(matids);
    fPressureFineMesh->ExpandSolution();
    
#ifdef PZDEBUG
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
    
    int64_t nc = cmeshPressure->NConnects();
    for (int64_t ic=nskeletonconnects; ic<nc; ic++) {
        cmeshPressure->ConnectVec()[ic].SetLagrangeMultiplier(1);
    }
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
            TPZMatLaplacian *matl2 = new TPZMatLaplacian((matid));
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
        TPZMatLaplacian *mathybrid = new TPZMatLaplacian(fPressureSkeletonMatId);
        mathybrid->SetDimension(fGMesh->Dimension()-1);
        fPressureFineMesh->InsertMaterialObject(mathybrid);
    }
    else
    {
        DebugStop();
    }

}


void TPZMHMixedMeshControl::CreateHDivPressureMHMMesh()
{
    TPZManVector<TPZCompMesh *,2 > cmeshes(2);
    cmeshes[0] = fFluxMesh.operator->();
    cmeshes[1] = fPressureFineMesh.operator->();
    TPZGeoMesh *gmesh = cmeshes[0]->Reference();
    if(!gmesh)
    {
        std::cout<< "Geometric mesh doesn't exist" << std::endl;
        DebugStop();
    }
    int dim = gmesh->Dimension();
    gmesh->ResetReference();
    // Malha computacional
    TPZCompMesh * MixedFluxPressureCmesh = fCMesh.operator->();
    
    
    
    MixedFluxPressureCmesh->SetDimModel(dim);
    MixedFluxPressureCmesh->SetAllCreateFunctionsMultiphysicElem();
    
    BuildMultiPhysicsMesh();
    TPZManVector<TPZCompMesh * ,2> meshvector(2);
    
    if(0)
    {
        std::ofstream out2("gmesh.txt");
        gmesh->Print(out2);
        std::ofstream out3("HDivMesh.txt");
        fFluxMesh->Print(out3);
        std::ofstream out4("PressureMesh.txt");
        fPressureFineMesh->Print(out4);
    }
    
    meshvector[0] = cmeshes[0];
    meshvector[1] = cmeshes[1];
    
    JoinSubdomains(meshvector, MixedFluxPressureCmesh);
    
    // Transferindo para a multifisica
    TPZBuildMultiphysicsMesh::AddElements(meshvector, MixedFluxPressureCmesh);
    TPZBuildMultiphysicsMesh::AddConnects(meshvector, MixedFluxPressureCmesh);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector, MixedFluxPressureCmesh);

#ifdef PZDEBUG
    if(1)
    {
        std::ofstream file("cmeshmphys.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(MixedFluxPressureCmesh, file,true);
        std::ofstream out("cmeshmphys.txt");
        MixedFluxPressureCmesh->Print(out);
    }
#endif
    std::pair<int,int> skelmatid(fSkeletonMatId,fSecondSkeletonMatId);
    CreateMultiPhysicsInterfaceElements(fGMesh->Dimension()-1);
    CreateMultiPhysicsInterfaceElements(fGMesh->Dimension()-2);
    MixedFluxPressureCmesh->CleanUpUnconnectedNodes();
    
    if(0)
    {
        std::ofstream out("multiphysics.txt");
        MixedFluxPressureCmesh->Print(out);
    }
    
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
void TPZMHMixedMeshControl::CreateInternalFluxElements()
{
    TPZGeoMesh *gmesh = fGMesh.operator->();
    int meshdim = gmesh->Dimension();
    TPZCompMesh * cmeshHDiv = fFluxMesh.operator->();
    gmesh->ResetReference();
    cmeshHDiv->LoadReferences();
    cmeshHDiv->SetDimModel(meshdim);
    cmeshHDiv->ApproxSpace().SetAllCreateFunctionsHDiv(meshdim);
    cmeshHDiv->SetDefaultOrder(fpOrderInternal);
    
    //Criar elementos computacionais malha MHM
    
    TPZGeoEl *gel = NULL;
    TPZGeoEl *gsubel = NULL;
    fConnectToSubDomainIdentifier[cmeshHDiv].Expand(10000);
    int64_t nel = fGMesh->NElements();
    
    for (auto it = fMHMtoSubCMesh.begin(); it != fMHMtoSubCMesh.end(); it++)
    {
        fGMesh->ResetReference();
        for (int64_t el = 0; el<nel; el++) {
            TPZGeoEl *gel = fGMesh->Element(el);
            if (!gel || gel->HasSubElement() || fGeoToMHMDomain[el] != it->first) {
                continue;
            }
            int64_t index;
            // create the flux element
            cmeshHDiv->CreateCompEl(gel, index);
            TPZCompEl *cel = cmeshHDiv->Element(index);
            /// associate the connects with the subdomain
            SetSubdomain(cel, it->first);
        }
        
        
    }

    fGMesh->ResetReference();
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
    std::map<int64_t, std::pair<int64_t,int64_t> >::iterator it = fInterfaces.begin();
    while (it != fInterfaces.end()) {
        int64_t elindex = it->first;
        // skip the boundary elements
        //        if (elindex == it->second.second) {
        //            it++;
        //            continue;
        //        }
        TPZGeoEl *gel = fGMesh->ElementVec()[elindex];
        int64_t index;
        // create an element to model the flux
        fFluxMesh->CreateCompEl(gel, index);
        TPZCompEl *cel = fFluxMesh->ElementVec()[index];
        int Side = gel->NSides()-1;
        TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
        SetSubdomain(cel, -1);
        
        if (elindex == it->second.second) {
            // this is a boundary element
            // set the side orientation of the boundary elements
            intel->SetSideOrient(Side, 1);
            SetSubdomain(cel, it->second.first);
        }
        else
        {
            if (it->second.first < it->second.second) {
                // set the flux orientation depending on the relative value of the element ids
                intel->SetSideOrient(Side, 1);
            }
            else
            {
                intel->SetSideOrient(Side, -1);
            }
            // this element will not be put in a subdomain
            SetSubdomain(cel, -1);
        }
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
        int64_t elindex = it->first;

        std::map<int64_t,std::list<TPZCompElSide> > subels;
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

    // for each pressure element create two interface elements
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
        bool keeplagrange = true;
        TPZCompMeshTools::CreatedCondensedElements(subcmesh, keeplagrange);
        subcmesh->CleanUpUnconnectedNodes();
        int numthreads = 0;
        int preconditioned = 0;
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
    TPZManVector<int64_t> shouldcreate(fGMesh->NElements(),0);
    std::set<int> matids;
    for (auto it : fCMesh->MaterialVec()) {
        matids.insert(it.first);
    }
    int64_t nel = fFluxMesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = fFluxMesh->Element(el);
        TPZGeoEl *gel = cel->Reference();
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
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = fGMesh->Element(el);
        if (shouldcreate[el])
        {
            int64_t index;
            fCMesh->CreateCompEl(gel, index);
        }
    }
}

