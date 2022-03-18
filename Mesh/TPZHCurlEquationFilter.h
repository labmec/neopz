//
// A filter to HCurl spaces.
//

#ifndef HCURL_EQUATION_FILTER_H
#define HCURL_EQUATION_FILTER_H

// TODO add doc

#include <iostream>
#include <map>
#include "pzstack.h"
#include <TPZVTKGeoMesh.h>

template <class T>
class TPZVec;

class TPZCompMesh;
class TPZGeoMesh;
class TPZCompElSide;
class TPZInterpolatedElement;
class TPZGeoEl;

template <class T, int N>
class TPZStack;

template <class TVar>
class TPZHCurlEquationFilter
{
       
public:
    // DATA STRUCTURES
    // 1 - VERTEX DATA STRUCTURE
    enum VertexStatusType {ENotTreatedVertex, ETreatedVertex};
    struct VertexFilter
    {
        // The edges connected to a vertex
        std::set<int> edge_connect;
        // The number of free edges, i.e, not removed in a vertex
        int64_t free_edges;
        // Edge associated removed
        int64_t removed_edge;
        // The vertex status: false - not treated; true - treated
        VertexStatusType status = ENotTreatedVertex;
    };

    // 2 - EDGE DATA STRUCTURE
    enum EdgeStatusType {EFreeEdge,ERemovedEdge,EBlockedEdge};
    struct EdgeFilter
    {
        //The edge index
        int64_t index;
        //The vertices connected to an edge
        std::set<int> vertex_connect;
        //The faces connected to an edge
        std::set<int> face_connect;
        //The edge status
        EdgeStatusType status = EFreeEdge;
        //Number of faces removed from the edge
        int64_t faces_removed = 0;
        //If it was removed, which is the vertex associated
        int64_t vertex_treated;
    };

    // 3 - FACE DATA STRUCTURE
    enum FaceStatusType {EFreeFace,ECryticalFace,EBlockedFace,EFortunateFace};
    struct FaceFilter
    {
        //The face index
        int64_t index;

        //The connects of a face
        std::set<int> edge_connect;

        //Face status:
        //  0 - free: 0 edges removed
        //  1 - crytical: 1 edge removed
        //  2 - blocked: 2 edges removed
        //  3 - fortunate: 1 edge blocked
        FaceStatusType status = EFreeFace;
    };

private:
    const int edgeDim{1};

    int fNumEdgesPerFace = 0;

    // 1 - Vertex Data Structure
    std::map<int64_t, VertexFilter> mVertex;
    // 2 - Edge Data Structure
    std::map<int64_t, EdgeFilter> mEdge;
    // 3 - Face Data Structure
    std::map<int64_t, FaceFilter> mFace;
    // 4 - The number of free and the corresponding treated nodes
    std::map<int64_t, std::set<int64_t>,std::greater<int>> freeEdgesToTreatedNodes;

    //The removed edges
    std::set<int64_t> removed_edges;

public:
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
    bool FilterEdgeEquations(TPZCompMesh* cmesh, TPZVec<int64_t> &activeEquations, bool domainHybridization = false);

    /**
        @brief Initialize the data structures, which are based only on the geometric mesh
        @param[in] gmesh geometric mesh.
    */
    void InitDataStructures(TPZGeoMesh *gmesh);

    /**
        @brief Associate an edge to a vertex and set it to be removed
        @param[out] treatVertex treated vertex
        @param[out] remEdge removed edge
    */
    void ChooseVertexAndEdge(int64_t &treatVertex, int64_t &remEdge);

    /**
        @brief Updates the edge and face status into the data structures
        @param[in] remEdge removed edge 
    */
    void UpdateEdgeAndFaceStatus(int64_t &remEdge);

    /**
        @brief Checks if all edges of a face have been removed. Returns a DebugStop() if true.
    */
    void CheckFaces();

    /**
        @brief Chooses the first edge to be removed: the first edge is the one
        that have the highest number of neighbour free edges
    */
    void FirstEdge();

    /**
        @brief Returns a set containing the indexes of all removed edges
    */
    std::set<int64_t> &GetRemovedEdges() {return removed_edges;};
    
    /**
        @brief Returns the vertex data structure
    */
    std::map<int64_t, VertexFilter> &GetVertexDataStructure() {return mVertex;};

    /**
        @brief Returns the edge data structure
    */
    std::map<int64_t, EdgeFilter> &GetEdgeDataStructure() {return mEdge;};
};
#endif