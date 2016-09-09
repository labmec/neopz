//
//  TPBRWellBBox.h
//  PZ
//
//  Created by Philippe Devloo on 3/27/15.
//
//

#ifndef __PZ__TPBRWellBBox__
#define __PZ__TPBRWellBBox__

#include <stdio.h>
#include <set>
#include "pzmanvector.h"
#include "tpzautopointer.h"
#include "pzgmesh.h"
#include "pzgeoelside.h"


/**
 * @brief Class which generates a mesh around the wellbore
 * @author Geo-Adapta
 */
class TPBRWellBBox
{
    /// Dimensions of the box
    REAL fLx, fLy, fZBottom, fZTop;
    
    /// Position of the well divisions
    TPZManVector<REAL> fWellPosY;
    
    /// Y Position of the outer nodes
    TPZManVector<REAL> fOuterPosY;
    
    /// Number of divisions towards the well
    int fNumDiv;
    
    /// Well diameter
    REAL fWellDiameter;
    
    /// Plane z coordinate (zero if the box is not cut
    REAL fPlaneZcoordinate;
    
    /// Geometric mesh generated
    TPZAutoPointer<TPZGeoMesh> fGMesh;
    
    struct TFace
    {
        TFace() : fOuterNodeIndices(4),fWellNodeIndices(4),fArcNodes(4)
        {
            
        }
        
        ~TFace()
        {
            
        }
        
        TFace(const TFace &copy)
        {
            DebugStop();
        }
        
        TFace &operator=(const TFace &copy)
        {
            DebugStop();
            return *this;
        }
        TPZManVector<int,4> fOuterNodeIndices;
        TPZManVector<int,4> fWellNodeIndices;
        
        TPZManVector<int, 4> fArcNodes;
        
    };
    
    /// indices of geometric nodes associated with each layer
    TPZManVector<TFace> fNodeLayers;
    
    /// add geometric nodes
    void GenerateNodeLayers();
    
    /// generate 4 reservoir elements and one well element between both faces
    void GenerateElementLayer(TFace &NodesLeft, TFace &NodesRight, bool createArc3DAtLeft, bool createArc3DAtRight);
    
    /// add to left and right elements
    void AddElementClosure();
    
    /// add the circular elements
    void GenerateArcElements();
    
    /// apply the successive refinements
    void ApplyDivisions();
    
    /// determine at which height the nodes are which need to be moved
    REAL HeightOfNearestPlane();
    
    /// move nodes at a certain height to a given height along their ribs
    void MovePoints(REAL fromheight, REAL toheight);
    
    /// divide elements whose edges cross the plane Z coordinate. return the set of newly created node indices
    void DirectionalDivide(std::list<TPZGeoElSide> &dividedsides);
    
    /// move the newly created nodes to the plane Z coordinate
    void MoveNodes(std::set<int> &newnodes);
    
    
    /// move the midsidenode of the element to the fPlaneZcoordinate
    void MoveMidSideNode(TPZGeoEl *gel, int side);
    
public:
    
    TPBRWellBBox(REAL WellDiam, REAL Lx, REAL Ly, REAL ZBottom, REAL ZTop) : fLx(Lx), fLy(Ly), fZBottom(ZBottom), fZTop(ZTop), fNumDiv(0),
    fWellDiameter(WellDiam), fPlaneZcoordinate(0.)
    {
        
    }
    
    TPBRWellBBox(const TPBRWellBBox&copy)
    {
        DebugStop();
    }
    
    TPBRWellBBox &operator=(const TPBRWellBBox &copy)
    {
        DebugStop();
        return *this;
    }
    
    /// Set the well partitions and the rib indices
    void SetWellDivision(TPZVec<REAL> &position, int divisiontowardswell);
    
    /// Set the coordinate of the plane which will cut the box
    void SetCutPlane(REAL planeZcoordinate)
    {
        fPlaneZcoordinate = planeZcoordinate;
    }
    
    /// cut the elements at the height of fPlaneZcoordinate
    void CutMeshatPlane();
    
    /// Generate the mesh based on the internal data structure
    TPZAutoPointer<TPZGeoMesh> GenerateMesh();
};

#endif /* defined(__PZ__TPBRWellBBox__) */
