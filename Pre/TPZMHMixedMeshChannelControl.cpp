//
//  TPZMHMixedMeshControl.cpp
//  PZ
//


#include "TPZMHMixedMeshChannelControl.h"
#include "TPZVecL2.h"
#include "pzbndcond.h"
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

#ifdef LOG4CXX
static PZLogger logger("pz.mhmixedmeshcontrol");
#endif


void TPZMHMixedMeshChannelControl::BuildComputationalMesh(bool usersubstructure,bool OpenChannel,std::map<int,std::pair<TPZGeoElSide,TPZGeoElSide>> oppen_channel){
 
    
    if (fpOrderInternal == 0 || fpOrderSkeleton == 0) {
        DebugStop();
    }
    InsertPeriferalMaterialObjects();
    CreateHDivMHMMesh();
    
    if(OpenChannel==true){
    for(auto it:oppen_channel){
        TPZCompEl *cel_from = it.second.first.Element()->Reference();
        int side_from = it.second.first.Side();
        TPZCompEl *cel_to = it.second.second.Element()->Reference();
        int side_to = it.second.second.Side();
        int c_from_index = cel_from->Index();
        int c_to_index = cel_to->Index();
        int connec_from_index = side_from - 4;
        int connec_to_index = side_to - 4;
        FluxMesh()->Element(c_from_index)->Connect(connec_from_index).RemoveDepend();
        FluxMesh()->Element(c_to_index)->Connect(connec_to_index).RemoveDepend();
        int conenec_from_to = FluxMesh()->Element(c_from_index)->ConnectIndex(connec_from_index);
        FluxMesh()->Element(c_to_index)->SetConnectIndex(connec_to_index, conenec_from_to);
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
