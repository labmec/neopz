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

        // int64_t UpdateFace(int removedEdge, std::set<int> faceEdges){

        //     return -1 = nao bloquead, senão connect index do edge bloqueado.
        // };
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

public:
    /**
         @brief Removes some equations associated with edges to ensure that
        the gradient of the lowest order H1 functions cannot be represented.
        @param[in] cmesh Computational mesh.
        @param[out] indices of all remaining equations.
        @return 0 if no errors were detected, 1 if a vertex was left untreated,
        2 if a vertex had all the adjacent edges removed.
    */
    bool FilterEdgeEquations(TPZCompMesh* cmesh, TPZVec<int64_t> &activeEquations, bool &domainHybridization, std::set<int64_t> &removed_edges);

    void InitDataStructures(TPZGeoMesh *gmesh);

    void ChooseNodeAndEdge(int64_t &treatNode, int64_t &remEdge);

    void UpdateEdgeAndFaceStatus(int64_t &remEdge);

    void CheckFaces();

    void FirstEdge();

    std::map<int64_t, VertexFilter> &GetVertexDataStructure() {return mVertex;};
    std::map<int64_t, EdgeFilter> &GetEdgeDataStructure() {return mEdge;};
};
#endif