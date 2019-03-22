//
//  TPZMultiphysicsCompMesh.cpp
//  pz
//
//  Created by Omar Durán on 3/21/19.
//

#include "TPZMultiphysicsCompMesh.h"
#include "pzmultiphysicscompel.h"

TPZMultiphysicsCompMesh::TPZMultiphysicsCompMesh() : TPZCompMesh(){
    
    m_active_approx_spaces.Resize(0);
    m_mesh_vector.Resize(0);
}

TPZMultiphysicsCompMesh::TPZMultiphysicsCompMesh(TPZGeoMesh * gmesh) : TPZCompMesh(gmesh){
    
    m_active_approx_spaces.Resize(0);
    m_mesh_vector.Resize(0);
}

TPZMultiphysicsCompMesh::TPZMultiphysicsCompMesh(const TPZMultiphysicsCompMesh &other) : TPZCompMesh(other) {
    m_active_approx_spaces  = other.m_active_approx_spaces;
    m_mesh_vector           = other.m_mesh_vector;
}

TPZMultiphysicsCompMesh & TPZMultiphysicsCompMesh::operator=(const TPZMultiphysicsCompMesh &other){
    
    if (this != & other) // prevent self-assignment
    {
        TPZCompMesh::operator=(other);
        m_active_approx_spaces  = other.m_active_approx_spaces;
        m_mesh_vector           = other.m_mesh_vector;
    }
    return *this;
}

void TPZMultiphysicsCompMesh::BuildMultiphysicsSpace(TPZManVector<int,5> & active_approx_spaces, TPZVec<TPZCompMesh * > & mesh_vector){
    
    if (m_mesh_vector.size() != m_active_approx_spaces.size()) {
        std::cout<< "TPZMultiphysicsCompMesh:: The vector provided should have the same size." << std::endl;
        DebugStop();
    }
    m_active_approx_spaces = active_approx_spaces;
    m_mesh_vector          = mesh_vector;
    int n_approx_spaces = m_mesh_vector.size();
    SetNMeshes(n_approx_spaces);
    
    SetAllCreateFunctionsMultiphysicElem();
    TPZCompMesh::AutoBuild();
    AdjustBoundaryElements();
    AddElements();
    AddConnects();
    LoadSolutionFromMeshes();
}

void TPZMultiphysicsCompMesh::AutoBuild(){
    
    DebugStop();
}

void TPZMultiphysicsCompMesh::AddElements(){

    TPZGeoMesh * geometry = Reference();
    geometry->ResetReference();
    int64_t n_cels = NElements();
    int n_approx_spaces = m_mesh_vector.size();
    for(int i_as = 0; i_as < n_approx_spaces; i_as++)
    {
        m_mesh_vector[i_as]->LoadReferences();
        int64_t icel;
        for(icel=0; icel < n_cels; icel++)
        {
            TPZCompEl * cel = ElementVec()[icel];
            TPZMultiphysicsElement * mfcel = dynamic_cast<TPZMultiphysicsElement *> (cel);
            TPZMultiphysicsInterfaceElement * mfint = dynamic_cast<TPZMultiphysicsInterfaceElement *>(cel);
            if(mfcel)
            {
                int64_t found = 0;
                TPZGeoEl * gel = mfcel->Reference();
                TPZStack<TPZCompElSide> celstack;
                TPZGeoElSide gelside(gel,gel->NSides()-1);
                
                if (gel->Reference()) {
                    mfcel->AddElement(gel->Reference(), i_as);
                    continue;
                }
                else
                {
                    TPZGeoEl *gelF = gel;
                    while(gelF->Father())
                    {
                        gelF = gelF->Father();
                        if (gelF->Reference()) {
#ifdef PZDEBUG
                            if (gelF->MaterialId() != gel->MaterialId()) {
                                DebugStop();
                            }
#endif
                            mfcel->AddElement(gelF->Reference(), i_as);
                            found = true;
                            break;
                        }
                    }
                }
                if (!found) {
                    mfcel->AddElement(0, i_as);
                }
            }
            else if (mfint) {
                //set up interface
            }
            else {
                DebugStop();
            }
        }
        geometry->ResetReference();
    }
    
    for (int64_t icel = 0; icel < n_cels; icel++) {
        TPZCompEl *cel = Element(icel);
        TPZMultiphysicsElement *mfel = dynamic_cast<TPZMultiphysicsElement *>(cel);
        if (!mfel) {
            continue;
        }
        mfel->SetActiveApproxSpaces(m_active_approx_spaces);
        mfel->InitializeIntegrationRule();
    }
    
}

void TPZMultiphysicsCompMesh::AddConnects(){
    
    int n_approx_spaces = m_mesh_vector.size();

    TPZVec<int64_t> FirstConnect(n_approx_spaces,0);
    int64_t nconnects = 0;
    for (int i_as = 0; i_as < n_approx_spaces; i_as++)
    {
        if (m_active_approx_spaces[i_as] == 0) {
            continue;
        }
        FirstConnect[i_as] = nconnects;
        nconnects += m_mesh_vector[i_as]->ConnectVec().NElements();
    }
    ConnectVec().Resize(nconnects);
    Block().SetNBlocks(nconnects);
    
    int64_t counter = 0;
    int64_t seqnum = 0;
    for (int i_as = 0; i_as < n_approx_spaces; i_as++)
    {
        if (m_active_approx_spaces[i_as] == 0) {
            continue;
        }
        
        int64_t ic;
        int64_t nc = m_mesh_vector[i_as]->ConnectVec().NElements();
        for (ic=0; ic<nc; ic++)
        {
            TPZConnect &refcon =  m_mesh_vector[i_as]->ConnectVec()[ic];
            ConnectVec()[counter] = refcon;
            if (refcon.SequenceNumber() >= 0) {
                ConnectVec()[counter].SetSequenceNumber(seqnum);
                ConnectVec()[counter].SetNState(refcon.NState());
                ConnectVec()[counter].SetNShape(refcon.NShape());
                ConnectVec()[counter].SetLagrangeMultiplier(refcon.LagrangeMultiplier());
                int ndof = refcon.NDof(*m_mesh_vector[i_as]);
                Block().Set(seqnum,ndof);
                seqnum++;
            }
            counter++;
        }

        for (ic=0; ic<nc; ic++)
        {
            TPZConnect &cn = ConnectVec()[FirstConnect[i_as]+ic];
            if (cn.HasDependency())
            {
                TPZConnect::TPZDepend *dep = cn.FirstDepend();
                while (dep) {
                    dep->fDepConnectIndex = dep->fDepConnectIndex+FirstConnect[i_as];
                    dep = dep->fNext;
                }
            }
        }
    }
    
    Block().SetNBlocks(seqnum);
    ExpandSolution();
    int64_t iel;
    int64_t nelem = NElements();
    for (iel = 0; iel < nelem; iel++)
    {
        
        
        TPZCompEl *celorig = ElementVec()[iel];
        TPZMultiphysicsElement *cel = dynamic_cast<TPZMultiphysicsElement *> (celorig);
        TPZMultiphysicsInterfaceElement *interfacee = dynamic_cast<TPZMultiphysicsInterfaceElement *> (celorig);
        if (interfacee) {
            continue;
        }
        if (!cel) {
            continue;
        }
        TPZStack<int64_t> connectindexes;
        int64_t imesh;
        std::list<TPZOneShapeRestraint> oneshape;
        for (int i_as = 0; i_as < n_approx_spaces; i_as++) {
            
            if (m_active_approx_spaces[i_as] == 0) {
                continue;
            }
            
            TPZCompEl *celref = cel->ReferredElement(i_as);
            if (!celref) {
                continue;
            }
            std::list<TPZOneShapeRestraint> celrest;
            celrest = celref->GetShapeRestraints();
            for (std::list<TPZOneShapeRestraint>::iterator it = celrest.begin(); it != celrest.end(); it++) {
                TPZOneShapeRestraint rest = *it;
                TPZOneShapeRestraint convertedrest(rest);
                for(int face = 0; face < rest.fFaces.size(); face++)
                {
                    int ic = rest.fFaces[face].first;
                    convertedrest.fFaces[face].first = ic+FirstConnect[i_as];
                }
                oneshape.push_back(convertedrest);
            }
            int64_t ncon = celref->NConnects();
            int64_t ic;
            for (ic=0; ic<ncon; ic++) {
                connectindexes.Push(celref->ConnectIndex(ic)+FirstConnect[i_as]);
            }
        }
        cel->SetConnectIndexes(connectindexes);
        for (std::list<TPZOneShapeRestraint>::iterator it = oneshape.begin(); it != oneshape.end(); it++) {
            cel->AddShapeRestraint(*it);
        }
    }
    
}

void TPZMultiphysicsCompMesh::LoadSolutionFromMeshes()
{
    int n_approx_spaces = m_active_approx_spaces.size();
    TPZManVector<int64_t> FirstConnectIndex(n_approx_spaces+1,0);
    for (int i_as = 0; i_as < n_approx_spaces; i_as++) {
        FirstConnectIndex[i_as+1] = FirstConnectIndex[i_as]+m_mesh_vector[i_as]->NConnects();
    }
    TPZBlock<STATE> &blockMF = Block();
    for (int i_as = 0; i_as < n_approx_spaces; i_as++) {
        
        if (m_active_approx_spaces[i_as] == 0) {
            continue;
        }
        
        int64_t ncon = m_mesh_vector[i_as]->NConnects();
        TPZBlock<STATE> &block = m_mesh_vector[i_as]->Block();
        int64_t ic;
        for (ic=0; ic<ncon; ic++) {
            TPZConnect &con = m_mesh_vector[i_as]->ConnectVec()[ic];
            int64_t seqnum = con.SequenceNumber();
            if(seqnum<0) continue;       /// Whether connect was deleted by previous refined process
            int blsize = block.Size(seqnum);
            TPZConnect &conMF = ConnectVec()[FirstConnectIndex[i_as]+ic];
            int64_t seqnumMF = conMF.SequenceNumber();
            for (int idf=0; idf<blsize; idf++) {
                blockMF.Put(seqnumMF, idf, 0, block.Get(seqnum, idf, 0));
            }
        }
    }
}

void TPZMultiphysicsCompMesh::LoadSolutionFromMultiPhysics()
{
    int n_approx_spaces = m_active_approx_spaces.size();
    TPZManVector<int64_t> FirstConnectIndex(n_approx_spaces+1,0);
    for (int i_as = 0; i_as < n_approx_spaces; i_as++) {
        if (m_active_approx_spaces[i_as] == 0) {
            continue;
        }
        FirstConnectIndex[i_as+1] = FirstConnectIndex[i_as]+m_mesh_vector[i_as]->NConnects();
    }
    TPZBlock<STATE> &blockMF = Block();
    for (int i_as = 0; i_as < n_approx_spaces; i_as++) {
        
        if (m_active_approx_spaces[i_as] == 0) {
            continue;
        }
        
        int64_t ncon = m_mesh_vector[i_as]->NConnects();
        TPZBlock<STATE> &block = m_mesh_vector[i_as]->Block();
        int64_t ic;
        for (ic=0; ic<ncon; ic++) {
            TPZConnect &con = m_mesh_vector[i_as]->ConnectVec()[ic];
            int64_t seqnum = con.SequenceNumber();
            if(seqnum<0) continue;       /// Whether connect was deleted by previous refined process
            int blsize = block.Size(seqnum);
            TPZConnect &conMF = ConnectVec()[FirstConnectIndex[i_as]+ic];
            int nelconnected = conMF.NElConnected();
            if (nelconnected == 0) {
                continue;
            }
            int64_t seqnumMF = conMF.SequenceNumber();
            int idf;
            for (idf=0; idf<blsize; idf++) {
                block.Put(seqnum, idf, 0, blockMF.Get(seqnumMF, idf, 0));
            }
        }
    }
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        sout << "Solutions of the referred meshes";
    }
#endif
    for (int i_as = 0; i_as < n_approx_spaces; i_as++) {
        if (m_active_approx_spaces[i_as] == 0) {
            continue;
        }
        m_mesh_vector[i_as]->LoadSolution(m_mesh_vector[i_as]->Solution());
    }
}
