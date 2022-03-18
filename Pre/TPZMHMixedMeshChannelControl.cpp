//
//  TPZMHMixedMeshControl.cpp
//  PZ
//


#include "TPZMHMixedMeshChannelControl.h"
#include "TPZBndCond.h"
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
#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.mhmixedmeshcontrol");
#endif


void TPZMHMixedMeshChannelControl::BuildComputationalMesh(bool usersubstructure,bool OpenChannel,std::map<int,std::pair<TPZGeoElSide,TPZGeoElSide>> oppen_channel){
 
    
    if (fpOrderInternal == 0 || fpOrderSkeleton == 0) {
        DebugStop();
    }
    InsertPeriferalMaterialObjects();
    CreateHDivMHMMesh();
    
    if(OpenChannel==true){
        for(auto it:oppen_channel){
            fGMesh->ResetReference();
            FluxMesh()->LoadReferences();
            TPZGeoEl *gel_from_orig = it.second.first.Element();
            int side_from = it.second.first.Side();
            TPZGeoEl *gel_to_orig = it.second.second.Element();
            int64_t gel_from_index = gel_from_orig->Index();
            int64_t gel_to_index = gel_to_orig->Index();
            int mhm1 = fGeoToMHMDomain[gel_from_index];
            int mhm2 = fGeoToMHMDomain[gel_to_index];
            if(mhm1 == mhm2) DebugStop();
            TPZGeoEl *gel_from = fGMesh->Element(gel_from_index);
            TPZGeoEl *gel_to = fGMesh->Element(gel_to_index);
            TPZCompEl *cel_from = gel_from->Reference();
            TPZCompEl *cel_to = gel_to->Reference();
            int side_to = it.second.second.Side();
            int c_from_index = cel_from->Index();
            int c_to_index = cel_to->Index();
            int connec_from_index = side_from - 4;
            int connec_to_index = side_to - 4;
            int64_t cfrom_index = cel_from->ConnectIndex(connec_from_index);
            int64_t cto_index = cel_to->ConnectIndex(connec_to_index);
            if(cfrom_index == cto_index) DebugStop();
            TPZConnect &cfrom = cel_from->Connect(connec_from_index);
            TPZConnect &cto = cel_to->Connect(connec_to_index);
            cfrom.RemoveDepend();
            cto.RemoveDepend();
    //        FluxMesh()->Element(c_from_index)->Connect(connec_from_index).RemoveDepend();
    //        FluxMesh()->Element(c_to_index)->Connect(connec_to_index).RemoveDepend();
            auto conenec_from_to = cel_from->ConnectIndex(connec_from_index);
            cel_to->SetConnectIndex(connec_to_index, conenec_from_to);
        }
    }
    

        
    
    fFluxMesh->CleanUpUnconnectedNodes();
    fFluxMesh->ExpandSolution();
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
        fCMesh->LoadReferences();
        std::ofstream out("cmesh_substruct.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(fCMesh->Reference(), out,true);
    }
#endif

}

/*
TPZMHMixedMeshControl::TPZMHMixedMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh, std::set<int64_t> &coarseindices) : TPZMHMeshControl(gmesh,coarseindices)
{
    fFluxMesh = new TPZCompMesh(gmesh);
    fFluxMesh->SetDimModel(gmesh->Dimension());
    fCMesh = new TPZCompMesh(gmesh);
    fCMesh->SetDimModel(gmesh->Dimension());
}
*/

void TPZMHMixedMeshChannelControl::HideTheElements()
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
    {
        int64_t ncon = fCMesh->NConnects();
        for (int64_t ic = 0; ic<ncon; ic++) {
            TPZConnect &c = fCMesh->ConnectVec()[ic];
            if(c.NElConnected() == 0 && c.HasDependency())
            {
//                std::cout << "Removing depend for connect " << ic << std::endl;
                c.RemoveDepend();
            }
        }
    }
    fCMesh->CleanUpUnconnectedNodes();
    
    GroupandCondenseElements();
    
    std::cout << "Finished substructuring\n";
}

int64_t TPZMHMixedMeshChannelControl::WhichSubdomain(TPZCompEl *cel)
{
    int ncon = cel->NConnects();
    std::set<int64_t> domains;
    TPZCompMesh *cmesh = cel->Mesh();
    TPZManVector<int64_t> &cvec = fConnectToSubDomainIdentifier[cmesh];
    for (int ic=0; ic<ncon; ic++)
    {
        int64_t cindex = cel->ConnectIndex(ic);
        if (cvec[cindex] != -1) {
            domains.insert(cvec[cindex]);
        }
    }
    //    if (domains.size() > 1) {
    //        for (int ic=0; ic<ncon; ic++) {
    //            int64_t cindex = cel->ConnectIndex(ic);
    //            std::cout << cindex << "|" << cvec[cindex] << " ";
    //        }
    //        std::cout << std::endl;
    //        DebugStop();
    //    }
    if (domains.size() ==0 || domains.size() ==2) {
        return -1;
    }
    int64_t domain = *domains.begin();
    return domain;
}
