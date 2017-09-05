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


// tototo
TPZMHMixedMeshControl::TPZMHMixedMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh, std::set<long> &coarseindices) : TPZMHMeshControl(gmesh,coarseindices)
{
    fFluxMesh = new TPZCompMesh(gmesh);
    fFluxMesh->SetDimModel(gmesh->Dimension());
    fCMesh = new TPZCompMesh(gmesh);
    fCMesh->SetDimModel(gmesh->Dimension());
}

TPZMHMixedMeshControl::TPZMHMixedMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh, TPZVec<long> &coarseindices) : TPZMHMeshControl(gmesh,coarseindices)
{
    fFluxMesh = new TPZCompMesh(gmesh);
    fFluxMesh->SetDimModel(gmesh->Dimension());
    fCMesh = new TPZCompMesh(gmesh);
    fCMesh->SetDimModel(gmesh->Dimension());
}

/// Insert Boundary condition objects that do not perform any actual computation
void TPZMHMixedMeshControl::InsertPeriferalMaterialObjects()
{
    TPZMaterial *mat = fCMesh->MaterialVec()[1];
    if (!mat) {
        DebugStop();
    }
    TPZFNMatrix<1,STATE> val1(1,1,0.), val2Flux(1,1,0.);
    int typePressure = 0;
    TPZBndCond * bcPressure = mat->CreateBC(mat, fPressureSkeletonMatId, typePressure, val1, val2Flux);
    //    bcN->SetForcingFunction(0,force);
    fCMesh->InsertMaterialObject(bcPressure);
    
    TPZBndCond * bcFlux = mat->CreateBC(mat, fSkeletonMatId, typePressure, val1, val2Flux);
    //    bcN->SetForcingFunction(0,force);
    fCMesh->InsertMaterialObject(bcFlux);
    
    TPZBndCond * bcSecondFlux = mat->CreateBC(mat, fSecondSkeletonMatId, typePressure, val1, val2Flux);
    //    bcN->SetForcingFunction(0,force);
    fCMesh->InsertMaterialObject(bcSecondFlux);
    
    int LagrangeMatIdLeft = 50;
    int LagrangeMatIdRight = 51;
    int nstate = 1;
    int dim = fGMesh->Dimension();
    TPZLagrangeMultiplier *matleft = new TPZLagrangeMultiplier(fLagrangeMatIdLeft,dim,nstate);
    TPZLagrangeMultiplier *matright = new TPZLagrangeMultiplier(fLagrangeMatIdRight,dim,nstate);
    
    fCMesh->InsertMaterialObject(matleft);
    fCMesh->InsertMaterialObject(matright);


}



/// Create all data structures for the computational mesh
void TPZMHMixedMeshControl::BuildComputationalMesh(bool usersubstructure)
{
    if (fpOrderInternal == 0 || fpOrderSkeleton == 0) {
        DebugStop();
    }
    CreateHDivMHMMesh();
    
    if (fHybridize) {
        Hybridize();
    }
    
    CreatePressureMHMMesh();
    
    
    CreateHDivPressureMHMMesh();
    std::cout << "Total number of equations " << fCMesh->Solution().Rows() << std::endl;
    fGlobalSystemWithLocalCondensationSize = fCMesh->NEquations();
    fGlobalSystemSize = fCMesh->Solution().Rows();

//    {
//        std::cout << "Connect/Subdomain ";
//        for (long i=0; i< fConnectToSubDomainIdentifier.size(); i++) {
//            std::cout << i << "|" << fConnectToSubDomainIdentifier[i] << "  ";
//            std::cout << std::endl;
//        }
//        std::ofstream out("multiphysics.txt");
//        fCMesh->Print(out);
//    }
    if (usersubstructure) {
        HideTheElements();
    }
    fNumeq = fCMesh->NEquations();

}


TPZCompMesh * TPZMHMixedMeshControl::CreateHDivMHMMesh()
{
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
    bc = matl2->CreateBC(matl2, -3, 0, val1, val2);
    cmeshHDiv->InsertMaterialObject(bc);
    bc = matl2->CreateBC(matl2, -4, 0, val1, val2);
    cmeshHDiv->InsertMaterialObject(bc);
    bc = matl2->CreateBC(matl2, fSkeletonMatId, 0, val1, val2);
    cmeshHDiv->InsertMaterialObject(bc);
    if(fSecondSkeletonMatId != 0)
    {
        bc = matl2->CreateBC(matl2, fSecondSkeletonMatId, 0, val1, val2);
        cmeshHDiv->InsertMaterialObject(bc);
    }
    if(fPressureSkeletonMatId != 0)
    {
        bc = matl2->CreateBC(matl2, fPressureSkeletonMatId, 0, val1, val2);
        cmeshHDiv->InsertMaterialObject(bc);
        TPZMatLaplacian *mathybrid = new TPZMatLaplacian(fPressureSkeletonMatId);
        mathybrid->SetDimension(gmesh->Dimension()-1);
        fPressureFineMesh->InsertMaterialObject(mathybrid);
    }
    // totototo
    bc = matl2->CreateBC(matl2, 105, 0, val1, val2);
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
    
    CreateInternalElements();
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
    return cmeshHDiv;
}

//void TPZMHMixedMeshControl::DuplicateNeighbouringConnects()
//{
//    TPZCompMesh *HDivMesh = fFluxMesh.operator->();
//    TPZGeoMesh *gmesh = HDivMesh->Reference();
//    int dimension = gmesh->Dimension();
//    gmesh->ResetReference();
//    HDivMesh->LoadReferences();
//    HDivMesh->ComputeNodElCon();
//    long nel = HDivMesh->NElements();
//    for (long el=0; el<nel; el++) {
//        TPZCompEl *cel = HDivMesh->Element(el);
//        TPZGeoEl *gel = cel->Reference();
//        if (!gel || gel->Dimension() != dimension) {
//            continue;
//        }
//        int nc = cel->NConnects();
//        for (int ic =0; ic<nc; ic++) {
//            TPZConnect &c = cel->Connect(ic);
//            if (c.HasDependency() && c.NElConnected() == 2)
//            {
//                // duplicate the connect
//                long cindex = HDivMesh->AllocateNewConnect(c);
//                TPZConnect &newc = HDivMesh->ConnectVec()[cindex];
//                newc = c;
//                c.DecrementElConnected();
//                newc.DecrementElConnected();
//                cel->SetConnectIndex(ic, cindex);
//            }
//        }
//    }
//    HDivMesh->ExpandSolution();
//}

void TPZMHMixedMeshControl::DuplicateNeighbouringConnects()
{
    TPZCompMesh *HDivMesh = fFluxMesh.operator->();
    TPZGeoMesh *gmesh = HDivMesh->Reference();
    int dimension = gmesh->Dimension();
    gmesh->ResetReference();
    HDivMesh->LoadReferences();
    HDivMesh->ComputeNodElCon();
    long nel = gmesh->NElements();
    for (long el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (!gel || gel->MaterialId() != fSkeletonMatId) {
            continue;
        }
        TPZCompEl *cel = gel->Reference();
        if(!cel)
        {
            continue;
        }
        TPZGeoElSide gelside(gel,gel->NSides()-1);
        TPZStack<TPZCompElSide> celstack;
        TPZGeoElSide neighbour = gelside.Neighbour();
        while (neighbour != gelside) {
            if (neighbour.Reference()) {
                celstack.Push(neighbour.Reference());
            }
            else
            {
                neighbour.HigherLevelCompElementList2(celstack, 0, 0);
            }
            neighbour = neighbour.Neighbour();
        }
        std::set<long> subdomains;
        for (long ist = 0; ist<celstack.size(); ist++) {
            long domain = WhichSubdomain(celstack[ist].Element());
            if (domain != -1)
            {
                subdomains.insert(domain);
            }
//            std::set<long> connects;
//            celstack[ist].Element()->BuildConnectList(connects);
//            for (std::set<long>::iterator it = connects.begin(); it != connects.end(); it++) {
//                if(fConnectToSubDomainIdentifier[*it] != -1 ) {
//                    subdomains.insert(fConnectToSubDomainIdentifier[*it]);
//                }
//            }
        }
        if (subdomains.size() > 2) {
            DebugStop();
        }
        if (subdomains.size() != 2) {
            continue;
        }
        for (long ist = 0; ist<celstack.size(); ist++) {
            int side = celstack[ist].Side();
            TPZCompEl *celst = celstack[ist].Element();
            std::set<long> connects;
            std::set<long> domains;
            long domain = WhichSubdomain(celst);
            if (domain == -1) {
                continue;
            }
//            celst->BuildConnectList(connects);
//            for (std::set<long>::iterator it = connects.begin(); it != connects.end(); it++) {
//                if(fConnectToSubDomainIdentifier[*it] != -1 ) {
//                    TPZConnect &c = HDivMesh->ConnectVec()[*it];
//                    if (!c.HasDependency()) {
//                        domains.insert(fConnectToSubDomainIdentifier[*it]);
//                    }
//                }
//            }
//            if (domains.size() == 0) {
//                continue;
//            }
//            if (domains.size() != 1) {
//                DebugStop();
//            }
//            long domain = *domains.begin();
            TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(celst);

            int nc = intel->NSideConnects(side);
            for (int ic =0; ic<nc; ic++) {
                TPZConnect &c = intel->SideConnect(ic,side);
                int sideconnectlocid = intel->SideConnectLocId(ic, side);
                long connectindex = intel->ConnectIndex(sideconnectlocid);
                if (c.HasDependency() && c.NElConnected() >1)
                {
                    // duplicate the connect
                    long cindex = HDivMesh->AllocateNewConnect(c);
                    TPZConnect &newc = HDivMesh->ConnectVec()[cindex];
                    newc = c;
                    c.DecrementElConnected();
                    newc.DecrementElConnected();
                    newc.SetSequenceNumber(cindex);
                    SetSubdomain(cindex, domain);
                    celst->SetConnectIndex(sideconnectlocid, cindex);
                }
                else if(c.NElConnected() > 1)
                {
                    // duplicate the connect
                    long newcindex = HDivMesh->AllocateNewConnect(c);
                    TPZConnect &newc = HDivMesh->ConnectVec()[newcindex];
                    newc = c;
                    c.DecrementElConnected();
                    SetSubdomain(newcindex, domain);
                    newc.DecrementElConnected();
                    newc.SetSequenceNumber(newcindex);
                    // put dependency
                    int nshape = c.NShape();
                    TPZFMatrix<REAL> dep(nshape,nshape);
                    dep.Identity();
                    newc.AddDependency(newcindex, connectindex, dep, 0, 0, nshape, nshape);
                    celst->SetConnectIndex(sideconnectlocid, newcindex);

                }
                else
                {
                    fConnectToSubDomainIdentifier[connectindex] = domain;
                }
            }
        }
    }
    HDivMesh->ExpandSolution();
    HDivMesh->CleanUpUnconnectedNodes();
}

TPZCompMesh * TPZMHMixedMeshControl::CreatePressureMHMMesh()
{
    TPZGeoMesh * gmesh = fGMesh.operator->();
    gmesh->ResetReference();
    
    long nskeletonconnects = fPressureFineMesh->NConnects();
    int porder = fpOrderInternal;
    TPZCompMesh * cmeshPressure = fPressureFineMesh.operator->();
    gmesh->ResetReference();
    cmeshPressure->SetName("PressureMesh");
    cmeshPressure->SetDimModel(gmesh->Dimension());
    cmeshPressure->ApproxSpace().SetAllCreateFunctionsContinuous();
    cmeshPressure->ApproxSpace().CreateDisconnectedElements(true);
    cmeshPressure->SetDefaultOrder(porder);
    TPZMatLaplacian *matl2 = new TPZMatLaplacian(1);
    matl2->SetDimension(3);
    cmeshPressure->InsertMaterialObject(matl2);
    std::set<int> matids;
    matids.insert(1);
    cmeshPressure->AutoBuild(matids);
    long nc = cmeshPressure->NConnects();
    long offset = fFluxMesh->NConnects();
    for (long ic=nskeletonconnects; ic<nc; ic++) {
        cmeshPressure->ConnectVec()[ic].SetLagrangeMultiplier(1);
    }
    gmesh->ResetReference();
    cmeshPressure->LoadReferences();
    std::map<long,long>::iterator it;
    for (it = fCoarseIndices.begin(); it != fCoarseIndices.end(); it++) {
        TPZGeoEl *gel = fGMesh->Element(it->first);
        TPZStack<TPZCompElSide> celstack;
        TPZGeoElSide gelside(gel,gel->NSides()-1);
        gelside.EqualorHigherCompElementList2(celstack, 0, 0);
        int nst = celstack.size();
        for (long ist=0; ist<nst; ist++) {
            SetSubdomain(celstack[ist].Element(), it->first,offset);
        }
    }
    
    return cmeshPressure;
}

TPZCompMesh * TPZMHMixedMeshControl::CreateHDivPressureMHMMesh()
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
    
    
    MixedFluxPressureCmesh->AutoBuild();
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
    
    // Transferindo para a multifisica
    TPZBuildMultiphysicsMesh::AddElements(meshvector, MixedFluxPressureCmesh);
    TPZBuildMultiphysicsMesh::AddConnects(meshvector, MixedFluxPressureCmesh);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector, MixedFluxPressureCmesh);

    if(0)
    {
        std::ofstream out("multiphysics.txt");
        MixedFluxPressureCmesh->Print(out);
    }

    if(0)
    {
        MixedFluxPressureCmesh->LoadReferences();
        
        // create the multiphysics interface elements
        long nel = fFluxMesh->NElements();
        for (long el=0; el<nel; el++) {
            TPZCompEl *cel = fFluxMesh->Element(el);
            if (!el) {
                continue;
            }
            TPZInterfaceElement *intface = dynamic_cast<TPZInterfaceElement *>(cel);
            if (!intface) {
                continue;
            }
            TPZGeoElSide gelleft,gelright;
            TPZCompElSide celleft,celright,celleftMF,celrightMF;
            celleft = intface->LeftElementSide();
            celright = intface->RightElementSide();
            gelleft = celleft.Reference();
            gelright = celright.Reference();
            celleftMF = gelleft.Reference();
            celrightMF = gelright.Reference();
            if(!celleftMF.Element() || ! celrightMF.Element())
            {
                DebugStop();
            }
            int matid = intface->Material()->Id();
            long index;
            new TPZMultiphysicsInterfaceElement(*MixedFluxPressureCmesh,intface->Reference(),index,celleftMF,celrightMF);
        }
    }

    CreateMultiPhysicsInterfaceElements(fFluxMesh->Dimension()-1);
    MixedFluxPressureCmesh->CleanUpUnconnectedNodes();
    
    return MixedFluxPressureCmesh;
    
}

void TPZMHMixedMeshControl::HideTheElements()
{
    bool KeepOneLagrangian = true;
    if (fHybridize) {
        KeepOneLagrangian = false;
    }
    typedef std::set<long> TCompIndexes;
    std::map<long, TCompIndexes> ElementGroups;
    TPZGeoMesh *gmesh = fCMesh->Reference();
    gmesh->ResetReference();
    int dim = gmesh->Dimension();
    fCMesh->LoadReferences();
    long nel = fCMesh->NElements();
    for (long el=0; el<nel; el++) {
        TPZCompEl *cel = fCMesh->Element(el);
        long domain = WhichSubdomain(cel);
        if (domain == -1) {
            continue;
        }
        ElementGroups[domain].insert(el);
    }
//    for (std::set<long>::iterator it= fCoarseIndices.begin(); it != fCoarseIndices.end(); it++) {
//        long iel = *it;
//        TPZGeoEl *gel = gmesh->Element(iel);
//        if (gel->Dimension() != dim && gel->MaterialId() > 0) {
//            DebugStop();
//        }
//        // we took any neighbour of gel and identified a mapindex with it??
//        TPZStack<TPZCompElSide> highlevel;
//        TPZGeoElSide gelside(gel,gel->NSides()-1);
//        gelside.HigherLevelCompElementList3(highlevel, 0, 0);
//        long nelst = highlevel.size();
//        for (long elst=0; elst<nelst; elst++) {
//            ElementGroups[iel].insert(highlevel[elst].Element()->Index());
//        }
//        if (gel->Reference()) {
//            if (nelst) {
//                DebugStop();
//            }
//            ElementGroups[iel].insert(gel->Reference()->Index());
//        }
//        int nsides = gel->NSides();
//        for (int is=0; is<nsides-1; is++) {
//            if (gel->SideDimension(is) != dim-1) {
//                continue;
//            }
//            TPZGeoElSide gelside(gel,is);
//            TPZStack<TPZCompElSide> highlevel;
//            gelside.EqualorHigherCompElementList3(highlevel, 0, 0);
//            long nelst = highlevel.size();
//            for (long elst=0; elst<nelst; elst++) {
//                if (highlevel[elst].Reference().Dimension() != dim-1 || highlevel[elst].Element()->Reference()->Dimension() != dim-1) {
//                    continue;
//                }
//                int matid = highlevel[elst].Element()->Reference()->MaterialId();
//                if(matid != fSkeletonMatId)
//                {
//                    ElementGroups[iel].insert(highlevel[elst].Element()->Index());
//                }
//            }
//
//        }
//    }
    if (ElementGroups.size() <= 10)
    {
        std::cout << "Number of element groups " << ElementGroups.size() << std::endl;
        std::map<long,TCompIndexes>::iterator it;
        for (it=ElementGroups.begin(); it != ElementGroups.end(); it++) {
            std::cout << "Group " << it->first << " group size " << it->second.size() << std::endl;
            std::cout << " elements ";
            std::set<long>::iterator its;
            for (its = it->second.begin(); its != it->second.end(); its++) {
                std::cout << *its << "|" << fCMesh->Element(*its)->Reference()->Index() << " ";
            }
            std::cout << std::endl;
        }
    }
    
    
    std::map<long,long> submeshindices;
    TPZCompMeshTools::PutinSubmeshes(fCMesh.operator->(), ElementGroups, submeshindices, KeepOneLagrangian);
    std::cout << "After putting in substructures\n";
    fCoarseIndices = submeshindices;
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
        for (std::map<long,long>::iterator it= fCoarseIndices.begin(); it != fCoarseIndices.end(); it++) {
            out << it->first << " " << it->second << " ";
        }
        out << std::endl;
    }
    /// indices of the skeleton elements and their left/right elements of the skeleton mesh
    out << "Skeleton element indices with associated left and right coarse element indices\n";
    {
        std::map<long, std::pair<long,long> >::iterator it;
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
    out << "Subdomain indices of the connects\n";
    for (long i=0; i<fConnectToSubDomainIdentifier.size(); i++) {
        if (i && !(i%20)) {
            out << std::endl;
        }
        out << "(" << i << "->" << fConnectToSubDomainIdentifier[i] << ") ";
    }
    out << std::endl;

    
}

void TPZMHMixedMeshControl::Hybridize()
{

    fGMesh->ResetReference();
    long nskel = fInterfaces.size();
    std::map<long, std::pair<long, long> >::iterator it;
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
        int side = gel->NSides()-1;
        
        TPZGeoElSide gelside(gel,side);
        TPZGeoElBC skeleton2(gelside,fSecondSkeletonMatId);
        fFluxMesh->SetAllCreateFunctionsHDiv();
        // create a flux boundary element
        long indexflux;
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
        TPZGeoElBC pressuregel(gelside,fPressureSkeletonMatId);
        // this will be changed to the pressure mesh
//        fFluxMesh->SetAllCreateFunctionsContinuous();
        fPressureFineMesh->SetAllCreateFunctionsContinuous();
        long indexpressure;
//        fFluxMesh->CreateCompEl(pressuregel.CreatedElement(), indexpressure);
        fPressureFineMesh->CreateCompEl(pressuregel.CreatedElement(), indexpressure);
        TPZInterpolatedElement *presel = dynamic_cast<TPZInterpolatedElement *>(fPressureFineMesh->Element(indexpressure));
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
    long count = 0;
    for (it = fInterfaces.begin(); it != fInterfaces.end(); it++) {
        TPZGeoEl *gel = fGMesh->Element(it->first);
        // skip the boundary elements
        if (it->first == it->second.second) {
            continue;
        }
        TPZCompEl *cel = gel->Reference();
        if (!cel || cel != fluxorig[count]) {
            DebugStop();
        }
        std::map<long,std::list<TPZCompElSide> > connectedmap;
        // connected elements will compute all the subelements of the skeleton element to left and right
        ConnectedElements(it->first, it->second, connectedmap);
        if (connectedmap.size() != 2) {
            DebugStop();
        }
        // identify the domains left and right of the skeleton element
        long origdepindex = fluxorig[count]->ConnectIndex(0);
        long newdepindex = fluxsecond[count]->ConnectIndex(0);
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
    long numflux = fluxorig.size();
    for (int iflux=0; iflux<numflux; iflux++) {
        long fluxgelindex = fluxorig[iflux]->Reference()->Index();
        int fluxmat = fluxorig[iflux]->Reference()->MaterialId();
        int pressmat = pressure[iflux]->Reference()->MaterialId();
        long pressuregelindex = pressure[iflux]->Reference()->Index();
        fluxorig[iflux]->SetReference(pressuregelindex);
        pressure[iflux]->SetReference(fluxgelindex);
        fluxorig[iflux]->Reference()->SetMaterialId(fluxmat);
        pressure[iflux]->Reference()->SetMaterialId(pressmat);
    }

    // create interface elements between the flux elements and the pressure element
    if(0)
    {
        fFluxMesh->LoadReferences();
        count = 0;
        for (it = fInterfaces.begin(); it != fInterfaces.end(); it++) {
            TPZGeoEl *gel = fGMesh->Element(it->first);
            // skip the boundary elements
            if (it->first == it->second.second) {
                continue;
            }
            TPZCompEl *cel = gel->Reference();
            if (!cel || cel != pressure[count]) {
                DebugStop();
            }
            int side = gel->NSides()-1;
            TPZCompElSide celpressureside(cel,side);
            TPZCompElSide fluxorigside(fluxorig[count],side);
            TPZCompElSide fluxsecondside(fluxsecond[count],side);
            TPZGeoElSide gelside(gel,side);
            TPZGeoElBC gbcleft(gelside,fLagrangeMatIdLeft);
            TPZGeoElBC gbcright(gelside,fLagrangeMatIdRight);
            long index1,index2;
            if (it->second.first < it->second.second) {
                new TPZInterfaceElement(fFluxMesh,gbcleft.CreatedElement(),index1,fluxorigside,celpressureside);
                new TPZInterfaceElement(fFluxMesh,gbcright.CreatedElement(),index2,fluxsecondside,celpressureside);
            }
            else
            {
                new TPZInterfaceElement(fFluxMesh,gbcright.CreatedElement(),index1,fluxorigside,celpressureside);
                new TPZInterfaceElement(fFluxMesh,gbcleft.CreatedElement(),index2,fluxsecondside,celpressureside);
            }
            count++;
        }
    }
    // The connects of the pressure mesh are not to be condensed
    fFluxMesh->ExpandSolution();
    fPressureFineMesh->ExpandSolution();
    long nc = fFluxMesh->NConnects();
    long ncpr = fPressureFineMesh->NConnects();
    for (long ic=0; ic<ncpr; ic++) {
        SetSubdomain(ic, -1, nc);
    }
    
}


// create the elements domain per domain with approximation spaces disconnected from each other
void TPZMHMixedMeshControl::CreateInternalElements()
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
    fConnectToSubDomainIdentifier.Expand(10000);
    
    // create the elements contained in the geometric domains associated with the coarse mesh
    for (std::map<long,long>::iterator it = fCoarseIndices.begin(); it != fCoarseIndices.end(); it++)
    {
        std::set<long> elset;
        bool LagrangeCreated = false;
        long iel = it->first;
        elset.insert(iel);
        //        std::map<long, std::pair<long,long> >::iterator it2;
        //        for (it2 = fInterfaces.begin(); it2 != fInterfaces.end(); it2++) {
        //            if (it2->first == it2->second.second && it2->second.first == *it) {
        //                elset.insert(it2->first);
        //            }
        //        }
        
        while (elset.size()) {
            std::set<long>::iterator itel = elset.begin();
            long elfirst = *itel;
            elset.erase(itel);
            gel = fGMesh->ElementVec()[elfirst];
            if(!gel) DebugStop();
            // we work with only the leaf elements
            if (! gel->HasSubElement())
            {
                long index;
                // create the flux element
                cmeshHDiv->CreateCompEl(gel, index);
                TPZCompEl *cel = cmeshHDiv->Element(index);
                /// associate the connects with the subdomain
                SetSubdomain(cel, it->first);
            }
            else
            {
                int nsubels = gel->NSubElements();
                for (int is=0; is<nsubels; is++)
                {
                    gsubel = gel->SubElement(is);
                    elset.insert(gsubel->Index());
                }
            }
        }
        
        fGMesh->ResetReference();
        
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
    std::map<long, std::pair<long,long> >::iterator it = fInterfaces.begin();
    while (it != fInterfaces.end()) {
        long elindex = it->first;
        // skip the boundary elements
        //        if (elindex == it->second.second) {
        //            it++;
        //            continue;
        //        }
        TPZGeoEl *gel = fGMesh->ElementVec()[elindex];
        long index;
        // create an element to model the flux
        fFluxMesh->CreateCompEl(gel, index);
        TPZCompEl *cel = fFluxMesh->ElementVec()[index];
        int Side = gel->NSides()-1;
        TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
        SetSubdomain(cel, -1);
        
        if (elindex == it->second.second) {
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
    // Apply restraints to the element/sides along the skeleton
    fFluxMesh->LoadReferences();
    it = fInterfaces.begin();
    while (it != fInterfaces.end()) {
        long elindex = it->first;

        std::map<long,std::list<TPZCompElSide> > subels;
        ConnectedElements(elindex, it->second, subels);

        TPZGeoEl *gel = fGMesh->ElementVec()[elindex];
        int side = gel->NSides()-1;
        TPZCompEl *cel = gel->Reference();
        // intel contains the flux skeleton element
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);

        for (std::map<long,std::list<TPZCompElSide> >::iterator itlist = subels.begin(); itlist != subels.end(); itlist++) {
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
std::map<long,long>::iterator it;
for(it = fCoarseIndices.begin(); it != fCoarseIndices.end(); it++)
{
    TPZGeoEl *gel = fGMesh->Element(it->first);
    TPZStack<TPZCompElSide> connected;
    TPZGeoElSide gelside(gel,gel->NSides()-1);
    gelside.EqualorHigherCompElementList2(connected, 0, 0);
    int nst = connected.size();
    for (long ist = 0; ist<nst; ist++) {
        TPZCompEl *cel = connected[ist].Element();
        SetSubdomain(cel, it->first);
    }
}
*/

/// Create the interfaces between the pressure elements of dimension dim
void TPZMHMixedMeshControl::CreateMultiPhysicsInterfaceElements(int dim)
{
    TPZCompMesh * cmeshPressure = fPressureFineMesh.operator->();
    
    // Computational multiplysics mesh
    TPZCompMesh * MixedFluxPressureCmesh = fCMesh.operator->();
    // load only dim dimensional elements
    MixedFluxPressureCmesh->Reference()->ResetReference();
    long nel = MixedFluxPressureCmesh->NElements();
    for (long el=0; el<nel; el++) {
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
    for (long el = 0; el<nel; el++) {
        TPZCompEl *cel = cmeshPressure->Element(el);
        if (!cel || !cel->Reference()) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (gel->Dimension() != dim) {
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
            if (neighbour.Element()->MaterialId() == fSkeletonMatId || neighbour.Element()->MaterialId() == fSecondSkeletonMatId) {
                celstack.Push(neighbour.Reference());
            }
            neighbour = neighbour.Neighbour();
        }
        if (celstack.size() != 2) {
            DebugStop();
        }
        TPZCompElSide celside = gelside.Reference();
        // create an interface between the multiphysics element associated with the pressure element
        // and the neighbouring elements
        TPZGeoElBC gbcleft(gelside, fLagrangeMatIdLeft);
        TPZGeoElBC gbcright(gelside, fLagrangeMatIdRight);
        long index;
//        celside.Element()->Print();
//        celstack[0].Element()->Print();
//        celstack[1].Element()->Print();
        new TPZMultiphysicsInterfaceElement(*MixedFluxPressureCmesh,gbcleft.CreatedElement(),index,celstack[0],celside);
        new TPZMultiphysicsInterfaceElement(*MixedFluxPressureCmesh,gbcright.CreatedElement(),index,celstack[1],celside);
    }
}

/// group and condense the elements
void TPZMHMixedMeshControl::GroupandCondenseElements()
{
    for (std::map<long,long>::iterator it=fCoarseIndices.begin(); it != fCoarseIndices.end(); it++) {
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

