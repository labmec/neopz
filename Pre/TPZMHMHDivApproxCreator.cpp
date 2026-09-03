#include "TPZMHMHDivApproxCreator.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzgeoelbc.h"
#include "TPZVTKGeoMesh.h"
#include "TPZNullMaterial.h"
#include "TPZNullMaterialCS.h"
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include "TPZAnalyticSolution.h"
#include "pzsubcmesh.h"
#include "pzcondensedcompel.h"
#include "TPZLagrangeMultiplierCS.h"
#include "pzconnect.h"
#include "TPZHDivApproxCreator.h"

TPZMHMHDivApproxCreator::TPZMHMHDivApproxCreator(TPZGeoMesh* gmesh) : TPZHDivApproxCreator(gmesh), TPZMHMApproxCreator() {
};


TPZMHMHDivApproxCreator::TPZMHMHDivApproxCreator(TPZGeoMesh* gmesh, TPZVec<int64_t>& elPartition) : TPZHDivApproxCreator(gmesh), TPZMHMApproxCreator(elPartition) {
};


TPZMultiphysicsCompMesh * TPZMHMHDivApproxCreator::BuildMultiphysicsCMesh(){

    this->CheckSetupConsistency();

    if (this->HybridType() != HybridizationType::ENone){
        this->ComputePeriferalMaterialIds();
        this->AddHybridizationGeoElements();
        std::ofstream out("GeoMeshHybrid.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(this->fGeoMesh, out);
    }

    int fNumMeshes = 4;
    TPZManVector<TPZCompMesh*,7> meshvec(fNumMeshes);
    int countMesh = 0;

    meshvec[countMesh++] = this->CreateHDivSpace();

    TPZNullMaterial<> *nullmat = new TPZNullMaterial<>(fSkeletonMatId,1,1);//Change to nstate when using Elasticity material
    meshvec[0]->InsertMaterialObject(nullmat);

    meshvec[0]->SetDefaultOrder(fPOrderSkeleton);
    std::set<int> matSkel={fSkeletonMatId};
    meshvec[0]->AutoBuild(matSkel);

    // skeleton elements have lagrange multiplier 4
    int64_t nel = meshvec[0]->NElements();
    for (int el = 0; el<nel; el++) {
        TPZCompEl *cel = meshvec[0]->Element(el);
        TPZGeoEl *gel = cel->Reference();
        if(gel->MaterialId() == fSkeletonMatId)
        {
            cel->Connect(0).SetLagrangeMultiplier(4);
        }
    }

    int lagLevelCounter = 1;
    meshvec[countMesh++] = this->CreateL2Space(this->GetDefaultOrder(),lagLevelCounter++);
    fDistFluxLevel = lagLevelCounter;
    meshvec[countMesh++] = this->CreateConstantSpace(lagLevelCounter++);
    fAvPresLevel = lagLevelCounter;
    meshvec[countMesh++] = this->CreateConstantSpace(lagLevelCounter++);    

    std::ofstream outflux("flux2.txt");
    meshvec[0]->Print(outflux);
    std::ofstream outpressure("pressure2.txt");
    meshvec[1]->Print(outpressure);
    std::ofstream outconstflux("constflux2.txt");
    meshvec[2]->Print(outconstflux);
    std::ofstream outconstpress("constpressure2.txt");
    meshvec[3]->Print(outconstpress);

    TPZMultiphysicsCompMesh *cmesh = this->CreateMultiphysicsSpace(meshvec);


    return cmesh;
}


void TPZMHMHDivApproxCreator::PutinSubstructures(TPZCompMesh &cmesh)
{
    // create subcompmeshes
    std::map<int64_t,TPZSubCompMesh *> submeshes;
    for(auto part : fElementPartition) {
        if(part == -1) continue;
        if(submeshes.find(part) == submeshes.end())
        {
            auto *sub = new TPZSubCompMesh(cmesh);
            submeshes[part] = sub;
        }
    }
    int64_t nel = cmesh.NElements();
    for (int64_t el = 0; el<nel; el++) {
        auto *cel = cmesh.Element(el);
        auto *sub = dynamic_cast<TPZSubCompMesh *>(cel);
        if(sub) continue;
        auto *gel = cel->Reference();
        auto index = gel->Index();
        auto part = fElementPartition[index];
        if(part >= 0) {
            auto *sub = submeshes[part];
            sub->TransferElement(&cmesh, el);
        }
    }
    cmesh.ComputeNodElCon();
    for(auto itsub : submeshes) {
        TPZCompEl *submesh = itsub.second;
        int64_t ncon = submesh->NConnects();
        for(int64_t ic = 0; ic<ncon; ic++) {
            auto &c = submesh->Connect(ic);
            if(c.LagrangeMultiplier() == fAvPresLevel) {
                c.IncrementElConnected();
                break;
            }
        }
    }
    for(auto itsub : submeshes) {
        TPZSubCompMesh *submesh = itsub.second;
        submesh->MakeAllInternal();
        submesh->CleanUpUnconnectedNodes();
    }
    cmesh.ComputeNodElCon();
    {
        int64_t ncon = cmesh.NConnects();
        for(int64_t ic = 0; ic<ncon; ic++) {
            auto &c = cmesh.ConnectVec()[ic];
            if(c.NElConnected() == 0 && c.HasDependency()) {
                c.RemoveDepend();
            }
        }
    }
    cmesh.CleanUpUnconnectedNodes();
    if(1)
    {
        std::ofstream out("cmesh_substruct.txt");
        cmesh.Print(out);
        std::ofstream out2("cmesh_substruct.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(&cmesh, out2);
    }
}

void TPZMHMHDivApproxCreator::CondenseElements(TPZCompMesh &cmesh)
{
    const bool real_sol = cmesh.GetSolType()!=ESolType::EComplex;
    cmesh.ComputeNodElCon();
    int64_t nel = cmesh.NElements();
    cmesh.ElementSolution().Resize(nel, 5);
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = cmesh.Element(el);
        auto *subcel = dynamic_cast<TPZSubCompMesh *>(cel);
        if(subcel) {
            CondenseElements(*subcel);
            subcel->CleanUpUnconnectedNodes();
            subcel->SetAnalysisSparse(0);
        //    subcel->SetAnalysisSkyline();
        }
        else if(cel)
        {
            TPZGeoEl *gel = cel->Reference();
            if(gel && gel->Dimension() == cmesh.Dimension()) {
                // prevent the pressure connect from being condensed
                int nc = cel->NConnects();
                for (int ic = 0; ic<nc; ic++) {
                    TPZConnect &c = cel->Connect(ic);
                    if(c.LagrangeMultiplier() == fAvPresLevel) {
                        c.IncrementElConnected();
                    }
                }
                if(real_sol){
                    new TPZCondensedCompElT<STATE>(cel);
                }else{
                    new TPZCondensedCompElT<CSTATE>(cel);
                }
            }
        }
    }
    auto *sub = dynamic_cast<TPZSubCompMesh *>(&cmesh);
    if(0 && !sub)
    {
        std::ofstream out("CondensedCompMesh.txt");
        cmesh.Print(out);
    }
}

