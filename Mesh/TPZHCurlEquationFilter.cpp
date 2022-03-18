//
// A filter to HCurl spaces.
//

#include "TPZHCurlEquationFilter.h"

#include "pzcmesh.h"
#include "pzgmesh.h"
#include <TPZVTKGeoMesh.h>

/**
    @brief Initialize the data structures, which are based only on the geometric mesh
    @param[in] gmesh geometric mesh.
*/
template <class TVar>
void TPZHCurlEquationFilter<TVar>::InitDataStructures(TPZGeoMesh* gmesh)
{
    
    //Loop over the elements
    for (auto gel : gmesh->ElementVec()){
    
        //Only tetrahedra is considered in the algorithm
        if (gel->Dimension() < 3) continue;
        
        const auto nEdges = gel->NSides(edgeDim);
        const auto firstEdge = gel->FirstSide(edgeDim);
        const auto firstFace = firstEdge + nEdges;
        const auto index = gel->Index();
        const auto ncorner = gel->NCornerNodes();
        const auto nsides = gel->NSides();
        if (gel->Type(nsides-1) == ETetraedro){
            fNumEdgesPerFace = 3;
        } else if (gel->Type(nsides-1) == ECube){
            fNumEdgesPerFace = 4;
        } else if (gel->Type(nsides-1) == EPrisma){
            fNumEdgesPerFace = 4;//Improve it
        }


        for(auto ie = firstEdge; ie < firstFace; ie++){
            // Vertex data structure            
            TPZGeoElSide edge_(gel,ie);
            const auto conIndex = edge_.Element()->Reference()->ConnectIndex(ie-firstEdge);
            
            //Store the connects of a given vertex 
            const int64_t v1 = edge_.SideNodeIndex(0);
            const int64_t v2 = edge_.SideNodeIndex(1);
            
            mVertex[v1].edge_connect.insert(conIndex);
            mVertex[v2].edge_connect.insert(conIndex);

            mEdge[conIndex].vertex_connect.insert(v1);
            mEdge[conIndex].vertex_connect.insert(v2);

            TPZStack<TPZGeoElSide> gelsides;
            gel->AllHigherDimensionSides(ie,2,gelsides);

            for (int i = 0; i < gelsides.size(); i++)
            {   
                // Edge data structure
                const auto con_face = edge_.Element()->Reference()->ConnectIndex(gelsides[i].Side()-ncorner);
                // Insert the face connect value.
                mEdge[conIndex].index = conIndex;           
                mEdge[conIndex].face_connect.insert(con_face);

                // Face data structure
                // Insert the edge connect value.
                mFace[con_face].index = conIndex;           
                mFace[con_face].edge_connect.insert(conIndex);
            }//i
        }//ie
    }//gel

    //Complete the data structures
    for(auto it = mVertex.begin(); it != mVertex.end(); ++it)
    // for (auto inode = 0; inode < mVertex.size(); inode++)
    {
        int size_edges = it->second.edge_connect.size();
        it->second.free_edges = size_edges;
    }

    //Cheks if an edge has more than 2 vertices
    for(auto it = mEdge.cbegin(); it != mEdge.cend(); ++it){
        if (it->second.vertex_connect.size() != 2){
            std::cout << "Edge with more than 2 vertices!\n";
            DebugStop();
        }
    }
    
    //Cheks if a face has more than 3 edges for Tetrahedron and 4 edges for Cube/Hexahedron
    // for(auto it = mFace.cbegin(); it != mFace.cend(); ++it){
    //     if (it->second.edge_connect.size() != fNumEdgesPerFace){
    //         std::cout << "Face with more than " << fNumEdgesPerFace << " edges!\n"; 
    //         DebugStop();
    //     }
    // }

}


/**
    @brief Associate an edge to a vertex and set it to be removed
    @param[out] treatVertex treated vertex
    @param[out] remEdge removed edge
*/
template <class TVar>
void TPZHCurlEquationFilter<TVar>::ChooseVertexAndEdge(int64_t &treatVertex, int64_t &remEdge)
{
    
    // The treated node must be connected to an already treated node by a free edge
    // Thus, we take a treated node and then a free edge and its 
    // corresponding second node to be the treated one
    auto highest = freeEdgesToTreatedNodes.rbegin()->first; // highest number of free edges
    int64_t base_node;
    int64_t ibase_free;
    bool flag = false;
    // We look at all the treated nodes, starting from the one with the highest number of free edges
    for (const auto& free_edges : freeEdgesToTreatedNodes){
        if (flag) break;
        // We loop over each treated vertex with the maximum value of free edges (called basis node) 
        for (auto ibase : free_edges.second){       
            if (flag) break;
            // We loop over each edge connected to the basis node
            for (auto iedge : mVertex[ibase].edge_connect){
                if (mEdge[iedge].status != EFreeEdge) continue; // only free edges can be removed
                if (flag) break;
                // For each edge connected to the basis node, check if the complementary 
                // node has been treated before (to avoid choosing the same treatedNode more than once)
                for (auto inode : mEdge[iedge].vertex_connect){
                    if (mVertex[inode].status == ETreatedVertex) {
                        //The complementary node has been tretated, go to the next one
                        continue;
                    }
                    // The node was not treated before, then it is the chosen one.
                    // Thus fill the variables with the proper values
                    base_node = ibase;
                    remEdge = iedge;
                    treatVertex = inode;
                    ibase_free = free_edges.first;
                    
                    flag = true;
                    break;
                }
            }
        }
    }

    mVertex[base_node].free_edges--;
    mVertex[treatVertex].free_edges--;
    mEdge[remEdge].status = ERemovedEdge;
    
    //Update the data structure.
    //Delete node and realocate the base node
    freeEdgesToTreatedNodes[ibase_free].erase(base_node);
    freeEdgesToTreatedNodes[ibase_free-1].insert(base_node);
    //Insert treatVertex
    freeEdgesToTreatedNodes[mVertex[treatVertex].free_edges].insert(treatVertex);

    //Clears the freeEdgesToTreatedNodes with 0 nodes in some free edges - to avoid a bug when looping std::map
    for (auto it = freeEdgesToTreatedNodes.cbegin(); it != freeEdgesToTreatedNodes.cend();) {
        if (it->second.size() == 0) {
            freeEdgesToTreatedNodes.erase(it++);
        } else {
            ++it;
        }
    }//it
}

/**
    @brief Updates the edge and face status into the data structures
    @param[in] remEdge removed edge 
*/
template <class TVar>
void TPZHCurlEquationFilter<TVar>::UpdateEdgeAndFaceStatus(int64_t &remEdge)
{
    std::set<int> blockedEdges;

    //Look at the faces where the removed edge is a part of and updates its status
    for (auto iface : mEdge[remEdge].face_connect)
    {
        auto status = mFace[iface].status;
        if (status == EFreeFace){
            mFace[iface].status = ECryticalFace;
        } else if (status == ECryticalFace || status == EFortunateFace) { // here there is a ?
            mFace[iface].status = EBlockedFace;

            //Loops over the edges of a face. If it is a free edge, then block it
            for (auto iedge : mFace[iface].edge_connect){
                if (mEdge[iedge].status == EFreeEdge){
                    mEdge[iedge].status = EBlockedEdge;
                    blockedEdges.insert(iedge);
                    break;
                }
            }
        } else if (status == EBlockedFace){
            DebugStop();
        } 
    }//iface

    // Now we can update all the faces in the mesh, as we can have some blocked 
    // edges and it will produce fortunate faces.
    for (auto iedge : blockedEdges)
    {
        for (auto iface : mEdge[iedge].face_connect){
            mFace[iface].status = EFortunateFace;
        }
    }    
}

/**
    @brief Checks if all edges of a face have been removed. Returns a DebugStop() if true.
*/
template <class TVar>
void TPZHCurlEquationFilter<TVar>::CheckFaces()
{

    //Check face
    for(auto it = mFace.cbegin(); it != mFace.cend(); ++it){
        int aux = 0; 
        for (auto iedge : it->second.edge_connect){
            if (mEdge[iedge].status == ERemovedEdge) aux++;
        }
        if (aux == fNumEdgesPerFace){
            std::cout << "Problem with face " << it->first << std::endl;
            DebugStop();
        }
    }//it
}

/**
    @brief Chooses the first edge to be removed: the first edge is the one
    that have the highest number of neighbour free edges
*/
template <class TVar>
void TPZHCurlEquationFilter<TVar>::FirstEdge()
{
    //start the algorithm by picking one edge
    //Here we chose the edge where the sum of free edges of both vertices is bigger
    int64_t firstEdge;
    int64_t maxFreeEdges = 0;
    for (auto iedge : mEdge)
    {
        int64_t localFreeEdges = 0;
        for (auto inode : iedge.second.vertex_connect)
        {
            localFreeEdges += mVertex[inode].free_edges;
        }
        if (localFreeEdges <= maxFreeEdges) continue;

        maxFreeEdges = localFreeEdges;
        firstEdge = iedge.first;
    }

    // Update the status of the chosen edge
    mEdge[firstEdge].status = ERemovedEdge;

    //Insert the chosen edge and the nodes to the freeEdgesToTreatedNodes structure;
    //In addition, mark the nodes as treated.
    for (auto inode : mEdge[firstEdge].vertex_connect)
    {
        mVertex[inode].free_edges--; //the node lose one free edge
        int64_t fEdges = mVertex[inode].free_edges;
        freeEdgesToTreatedNodes[fEdges].insert(inode);
        mVertex[inode].status = ETreatedVertex;

        mVertex[inode].removed_edge = firstEdge;
        mEdge[firstEdge].vertex_treated = inode;
    }
}


/**
    @brief Removes some equations associated with edges to ensure that
    the gradient of the lowest order H1 functions cannot be represented.
    @param[in] cmesh Computational mesh.
    @param[out] activeEquations a vector containing the index of all remaining equations.
    @param[in] domainHybridization If some kind of hybridization is applied to fluxes, this input 
    should be set as true as in this case the edge connects are disconnected and edges need to be 
    removed elementwise.
    @return 0 if no errors were detected, 1 if a vertex was left untreated,
    2 if a vertex had all the adjacent edges removed.
*/
template <class TVar>
bool TPZHCurlEquationFilter<TVar>::FilterEdgeEquations(TPZCompMesh* cmesh,
                    TPZVec<int64_t> &activeEquations, bool domainHybridization)
{

    cmesh->LoadReferences();
    const auto gmesh = cmesh->Reference();  
    auto nnodes = gmesh->NNodes();

    TPZVec<bool> done_vertices(nnodes, false);

    if (domainHybridization)
    {
        done_vertices = true;
        for (auto cel : cmesh->ElementVec())
        {   
            if (!cel) continue;
            cel->LoadElementReference();
            auto gel = cel->Reference();
            if (!gel) continue;
            int dimg = gel->Dimension();
            if (gel->Dimension() != 3) continue;

            switch (gel->Type())
            {
            case ETetraedro:
                removed_edges.insert(cel->ConnectIndex(0));
                removed_edges.insert(cel->ConnectIndex(1));
                removed_edges.insert(cel->ConnectIndex(5));
                break;
            case ECube:
                removed_edges.insert(cel->ConnectIndex(0));
                removed_edges.insert(cel->ConnectIndex(1));
                removed_edges.insert(cel->ConnectIndex(3));
                removed_edges.insert(cel->ConnectIndex(6));
                removed_edges.insert(cel->ConnectIndex(7));
                removed_edges.insert(cel->ConnectIndex(8));
                removed_edges.insert(cel->ConnectIndex(11));
                break;
            
            default:
                DebugStop();
                break;
            }
        }      
   
    } else {
        //Initialize the data structures
        InitDataStructures(gmesh);

        // Choose the first edge to be removed
        FirstEdge();

        // Always starts with 1 nodes treated + 1 node which edge will not be removed
        int64_t treated_nodes = 2;

        // While there is a node to be treated
        while (treated_nodes < nnodes)
        {
            // - Escolher um nó do mapa
            int64_t treatVertex, remEdge;
            ChooseVertexAndEdge(treatVertex,remEdge);

            mVertex[treatVertex].removed_edge = remEdge;
            mVertex[treatVertex].status = ETreatedVertex;
            mEdge[remEdge].vertex_treated = treatVertex;

            // Updates the edges (if some needs to be blocked) and faces status.
            // UpdateEdgeAndFaceStatus(remEdge);
            
            // Checks if all edges of a face have been removed
            CheckFaces();

            treated_nodes++;
            
        }//while

        for (size_t i = 0; i < nnodes; i++)//skip the first node since it is not associated to a removed edge
        {
            if (mVertex[i].status == ETreatedVertex) {
                done_vertices[i] = true;
                removed_edges.insert(mVertex[i].removed_edge);
            }
        }
    }
    
    // std::cout << "Removed (" << removed_edges.size() << ") = ";
    // for (auto rem : removed_edges)
    // {
    //     std::cout << rem << " " ;
    // }
    // std::cout << std::endl;

    // for (int i = 0; i < nnodes; i++)
    // {
    //     std::cout << "MVertex (" << i << ") = " << mVertex[i].removed_edge << std::endl;
    // }
    
    



    if (!domainHybridization){
        bool check_all_vertices{true};
        for(auto iv = 0; iv < nnodes; iv++){
            bool v = true;
            if (mVertex[iv].status != ETreatedVertex) v = false;
            check_all_vertices = check_all_vertices && v;
            if(!v){
            break;
            }
        }
        if (!check_all_vertices) DebugStop();
     
        bool check_edges_left{true};
        for(auto iv = 0; iv < nnodes; iv++){
            const auto all_edges = mVertex[iv].edge_connect;
            bool local_check{false};
            for(auto edge: all_edges){
                if(removed_edges.find(edge) == removed_edges.end()){
                    local_check = true;
                    break;
                }
            }
            check_edges_left = local_check && check_edges_left;
        }

        // if (!check_edges_left) DebugStop();
    }


    activeEquations.Resize(0);
    for (int iCon = 0; iCon < cmesh->NConnects(); iCon++) {
        if (removed_edges.find(iCon) == removed_edges.end()) {
            auto &con = cmesh->ConnectVec()[iCon];
            if (con.HasDependency()){
                continue;
            }
            const auto seqnum = con.SequenceNumber();
            if (seqnum<0) continue;
            if (con.IsCondensed()) continue;


            const auto pos = cmesh->Block().Position(seqnum);
            const auto blocksize = cmesh->Block().Size(seqnum);
            if (blocksize == 0){
                continue;
            }
            
            const auto vs = activeEquations.size();
            activeEquations.Resize(vs + blocksize);
            for (auto ieq = 0; ieq < blocksize; ieq++) {
                activeEquations[vs + ieq] = pos + ieq;
            }
        }
    }
    
    if (!domainHybridization){
        if (removed_edges.size() != nnodes-1){
            std::cout << "Removed " << removed_edges.size() << "/" << nnodes-1 << " allowed connects."  << std::endl;
            DebugStop();
        }
    } 

    return 0;
}


template class TPZHCurlEquationFilter<STATE>;
template class TPZHCurlEquationFilter<CSTATE>;


