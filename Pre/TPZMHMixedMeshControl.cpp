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
#include "mixedpoisson.h"

#include "pzsubcmesh.h"

#include "pzbuildmultiphysicsmesh.h"
#include "TPZCompMeshTools.h"

TPZMHMixedMeshControl::TPZMHMixedMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh, std::set<long> &coarseindices) : TPZMHMeshControl(gmesh,coarseindices)
{
    fFluxMesh = new TPZCompMesh(gmesh);
    fFluxMesh->SetDimModel(gmesh->Dimension());
}



/// Create all data structures for the computational mesh
void TPZMHMixedMeshControl::BuildComputationalMesh(bool usersubstructure)
{
    fFluxMesh = CreateHDivMHMMesh(fGMesh.operator->(),fpOrderInternal);
    DuplicateNeighbouringConnects(fFluxMesh.operator->());
    fPressureFineMesh = CreatePressureMHMMesh(fGMesh.operator->(),fpOrderInternal);
}

/// will create 1D elements on the interfaces between the coarse element indices
void TPZMHMixedMeshControl::CreateCoarseInterfaces(int matid)
{
    
}

TPZCompMesh * TPZMHMixedMeshControl::CreateHDivMHMMesh(TPZGeoMesh * gmesh, int porder)
{
    int meshdim = gmesh->Dimension();
    TPZCompMesh * cmeshHDiv = new TPZCompMesh(gmesh);
    cmeshHDiv->SetDimModel(meshdim);
    cmeshHDiv->ApproxSpace().SetAllCreateFunctionsHDiv(meshdim);
    cmeshHDiv->SetDefaultOrder(porder);
    TPZVecL2 *matl2 = new TPZVecL2(1);
    cmeshHDiv->InsertMaterialObject(matl2);
    matl2 = new TPZVecL2(2);
    cmeshHDiv->InsertMaterialObject(matl2);
    TPZFNMatrix<1,STATE> val1(1,1,0.),val2(1,1,0.);
    TPZBndCond *bc = matl2->CreateBC(matl2, -1, 0, val1, val2);
    cmeshHDiv->InsertMaterialObject(bc);
    bc = matl2->CreateBC(matl2, -2, 0, val1, val2);
    cmeshHDiv->InsertMaterialObject(bc);
    cmeshHDiv->AutoBuild();
    
#ifdef PZDEBUG
    {
        std::ofstream outmesh("BigHDivMesh.txt");
        cmeshHDiv->Print(outmesh);
    }
#endif
    return cmeshHDiv;
}

void TPZMHMixedMeshControl::DuplicateNeighbouringConnects(TPZCompMesh * HDivMesh)
{
    TPZGeoMesh *gmesh = HDivMesh->Reference();
    int dimension = gmesh->Dimension();
    gmesh->ResetReference();
    HDivMesh->LoadReferences();
    HDivMesh->ComputeNodElCon();
    long nel = HDivMesh->NElements();
    for (long el=0; el<nel; el++) {
        TPZCompEl *cel = HDivMesh->Element(el);
        TPZGeoEl *gel = cel->Reference();
        if (!gel || gel->Dimension() != dimension) {
            continue;
        }
        int nc = cel->NConnects();
        for (int ic =0; ic<nc; ic++) {
            TPZConnect &c = cel->Connect(ic);
            if (c.HasDependency() && c.NElConnected() == 2)
            {
                // duplicate the connect
                long cindex = HDivMesh->AllocateNewConnect(c);
                TPZConnect &newc = HDivMesh->ConnectVec()[cindex];
                newc = c;
                c.DecrementElConnected();
                newc.DecrementElConnected();
                cel->SetConnectIndex(ic, cindex);
            }
        }
    }
    HDivMesh->ExpandSolution();
}

TPZCompMesh * TPZMHMixedMeshControl::CreatePressureMHMMesh(TPZGeoMesh * gmesh, int porder)
{
    TPZCompMesh * cmeshPressure = new TPZCompMesh(gmesh);
    cmeshPressure->SetDimModel(gmesh->Dimension());
    cmeshPressure->ApproxSpace().SetAllCreateFunctionsContinuous();
    cmeshPressure->ApproxSpace().CreateDisconnectedElements(true);
    cmeshPressure->SetDefaultOrder(porder);
    TPZMatLaplacian *matl2 = new TPZMatLaplacian(1);
    matl2->SetDimension(3);
    cmeshPressure->InsertMaterialObject(matl2);
    cmeshPressure->AutoBuild();
    long nc = cmeshPressure->NConnects();
    for (long ic=0; ic<nc; ic++) {
        cmeshPressure->ConnectVec()[ic].SetLagrangeMultiplier(1);
    }
    return cmeshPressure;
}

TPZCompMesh * TPZMHMixedMeshControl::CreateHDivPressureMHMMesh(TPZVec<TPZCompMesh * > & cmeshes)
{
    TPZGeoMesh *gmesh = cmeshes[0]->Reference();
    if(!gmesh)
    {
        std::cout<< "Geometric mesh doesn't exist" << std::endl;
        DebugStop();
    }
    int dim = gmesh->Dimension();
    
    const int typeFlux = 1, typePressure = 0;
    TPZFMatrix<STATE> val1(1,1,0.), val2Flux(1,1,0.), val2Pressure(1,1,10.);
    val2Pressure(0,0) = 1000.;
    
    // Malha computacional
    TPZCompMesh * MixedFluxPressureCmesh = new TPZCompMesh(gmesh);
    
    // Material medio poroso
    TPZMixedPoisson * mat = new TPZMixedPoisson(1,dim);
    mat->SetSymmetric();
    //    mat->SetForcingFunction(One);
    MixedFluxPressureCmesh->InsertMaterialObject(mat);
    
    
    
    // Bc N
    TPZBndCond * bcN = mat->CreateBC(mat, -1, typePressure, val1, val2Pressure);
//    TPZAutoPointer<TPZFunction<STATE> > force = new TPZDummyFunction<STATE>(DirichletValidacao);
//    bcN->SetForcingFunction(0,force);
    MixedFluxPressureCmesh->InsertMaterialObject(bcN);
    
    // Bc S
    TPZBndCond * bcS = mat->CreateBC(mat, -2, typePressure, val1, val2Pressure);
//    bcS->SetForcingFunction(0, force);
    MixedFluxPressureCmesh->InsertMaterialObject(bcS);
    
    
    
    
    
    MixedFluxPressureCmesh->SetDimModel(dim);
    MixedFluxPressureCmesh->SetAllCreateFunctionsMultiphysicElem();
    MixedFluxPressureCmesh->AutoBuild();
    
    TPZManVector<TPZCompMesh * ,2> meshvector(2);
    
    
    meshvector[0] = cmeshes[0];
    meshvector[1] = cmeshes[1];
    
    // Transferindo para a multifisica
    TPZBuildMultiphysicsMesh::AddElements(meshvector, MixedFluxPressureCmesh);
    TPZBuildMultiphysicsMesh::AddConnects(meshvector, MixedFluxPressureCmesh);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector, MixedFluxPressureCmesh);
    
    return MixedFluxPressureCmesh;
    
}

void HideTheElements(TPZCompMesh * Multiphysics, bool KeepOneLagrangian, TPZVec<long> &coarseindices)
{
    typedef std::set<long> TCompIndexes;
    std::map<long, TCompIndexes> ElementGroups;
    TPZGeoMesh *gmesh = Multiphysics->Reference();
    gmesh->ResetReference();
    int dim = gmesh->Dimension();
    Multiphysics->LoadReferences();
    long nelg = coarseindices.NElements();
    for (long iel=0; iel<nelg; iel++) {
        long el = coarseindices[iel];
        TPZGeoEl *gel = gmesh->Element(el);
        if (gel->Dimension() != dim && gel->MaterialId() > 0) {
            DebugStop();
        }
        // we took any neighbour of gel and identified a mapindex with it??
        TPZStack<TPZCompElSide> highlevel;
        TPZGeoElSide gelside(gel,gel->NSides()-1);
        gelside.HigherLevelCompElementList3(highlevel, 0, 0);
        long nelst = highlevel.size();
        for (long elst=0; elst<nelst; elst++) {
            ElementGroups[el].insert(highlevel[elst].Element()->Index());
        }
        if (gel->Reference()) {
            if (nelst) {
                DebugStop();
            }
            ElementGroups[el].insert(gel->Reference()->Index());
        }
    }
    std::cout << "Number of element groups " << ElementGroups.size() << std::endl;
    std::map<long,TCompIndexes>::iterator it;
    for (it=ElementGroups.begin(); it != ElementGroups.end(); it++) {
        std::cout << "Group " << it->first << " group size " << it->second.size() << std::endl;
        std::cout << " elements ";
        std::set<long>::iterator its;
        for (its = it->second.begin(); its != it->second.end(); its++) {
            std::cout << *its << " ";
        }
        std::cout << std::endl;
    }
    
    std::set<long> submeshindices;
    TPZCompMeshTools::PutinSubmeshes(Multiphysics, ElementGroups, submeshindices, KeepOneLagrangian);
    std::cout << "After putting in substructures\n";
    
    Multiphysics->ComputeNodElCon();
    Multiphysics->CleanUpUnconnectedNodes();
    for (std::set<long>::iterator it=submeshindices.begin(); it != submeshindices.end(); it++) {
        TPZCompEl *cel = Multiphysics->Element(*it);
        TPZSubCompMesh *subcmesh = dynamic_cast<TPZSubCompMesh *>(cel);
        if (!subcmesh) {
            DebugStop();
        }
        TPZCompMeshTools::GroupElements(subcmesh);
        subcmesh->ComputeNodElCon();
        TPZCompMeshTools::CreatedCondensedElements(subcmesh, KeepOneLagrangian);
        subcmesh->CleanUpUnconnectedNodes();
        subcmesh->SetAnalysisSkyline(16, 0, 0);
    }
    //    Multiphysics->ComputeNodElCon();
    //    Multiphysics->CleanUpUnconnectedNodes();
    std::cout << "Finished substructuring\n";
}

