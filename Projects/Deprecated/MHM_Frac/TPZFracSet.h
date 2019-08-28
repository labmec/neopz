//
//  FracPre.cpp
//  PZ
//
//  Created by Philippe Devloo on 07/09/17.
//
//
#ifndef TPZFRACSETH
#define TPZFRACSETH

#include <stdio.h>
#include <list>
#include <fstream>
#include "pzmanvector.h"
#include "pzadmchunk.h"
#include "pzgnode.h"
#include "pzgmesh.h"



struct TPZFracture
{
    int fOrigId;
    
    int fMatId = -1;
    
    TPZManVector<int64_t,2> fNodes;
    
    REAL fFracPerm = 500;
    
    std::string fPhysicalName = "FRAC";
    
    TPZFracture() : fOrigId(-1)
    {
        
    }
    
    TPZFracture(int id, int matid, int64_t left, int64_t right) : fOrigId(id), fMatId(matid), fNodes(2), fFracPerm(100.)
    {
        fNodes[0] = left;
        fNodes[1] = right;
    }
    
    TPZFracture(const TPZFracture &copy) : fOrigId(copy.fOrigId), fMatId(copy.fMatId), fNodes(copy.fNodes), fFracPerm(copy.fFracPerm), fPhysicalName(copy.fPhysicalName)
    {
        
    }
    
    TPZFracture &operator=(const TPZFracture &copy)
    {
        fOrigId = copy.fOrigId;
        fMatId = copy.fMatId;
        fNodes = copy.fNodes;
        fFracPerm = copy.fFracPerm;
        fPhysicalName = copy.fPhysicalName;
        return *this;
    }
};

struct XOrder
{
    bool operator()(const TPZGeoNode &a, const TPZGeoNode &b) const
    {
        return a.Coord(0) < b.Coord(0);
    }
};

struct YOrder
{
    bool operator()(const TPZGeoNode &a, const TPZGeoNode &b) const
    {
        return a.Coord(1) < b.Coord(1);
    }
};
struct TPZFracSet
{
    
    const int matid_internal_frac = 1;
    
    const int matid_MHM_line = 50;
    
    const int matid_MHM_frac = 55;
    
    const int matid_BC = -1;
    
    TPZAdmChunkVector<TPZFracture> fFractureVec;

    TPZAdmChunkVector<TPZGeoNode> fNodeVec;
    
    TPZManVector<REAL> fMeshSizeAtNodes;
    
    std::map<uint64_t,int64_t> fPointMap;
    
    // ordered set of nodes which are along horizontal MHM lines (computed in AddMHMNodes)
    TPZManVector<std::map<REAL, int64_t> > fHorizontalLines;
    
    // ordered set of nodes which are along vertical MHM lines (computed in AddMHMNodes)
    TPZManVector<std::map<REAL, int64_t> > fVerticalLines;
    
    // set of node indices which are on MHM boundaries (computed in ADDMHMNodes)
    std::set<int64_t> fMHMNodes;
    
    /// NeoPZ representation of the fracture network
    TPZGeoMesh fgmesh;
    
    /// Domain size
    TPZManVector<REAL,3> fLowLeft, fTopRight;
    
    /// Number of MHM domains in x and y
    int fNumMHMDomains;
    /// Definition of MHM spacing in float values
    TPZManVector<REAL,3> fMHMSpacing;
    
    /// Definition of MHM spacing in integer values (multiple of fTol
    TPZManVector<uint64_t,3> fMHMSpacingInt;
    
    /// Distance tolerance
    REAL fTol = 1;
    
    /// Desired finite element size
    REAL fElementSize = 30.;
    
    /// Desired element size at tips
    REAL fMinElementSize = 1.;
    
    /// Name for the 2 dimensional physical name
    std::string fPhysicalname = "Darcy";
    
    TPZFracSet() : fLowLeft(3,0.), fTopRight(3,1000.), fMHMSpacing(2,31), fMHMSpacingInt(2,0), fTol(1.)
    {
        fNumMHMDomains = 31;
        SetMHMSpacing(fNumMHMDomains);
        SetTol(fTol);
    }
    
    TPZFracSet(const TPZFracSet &copy)
    {
        DebugStop();
    }
    
    TPZFracSet &operator=(const TPZFracSet &copy)
    {
        DebugStop();
        static TPZFracSet dummy;
        return dummy;
    }
    
    void SetMHMSpacing(int numdomains)
    {
        fNumMHMDomains = numdomains;
        fMHMSpacing[0] = (fTopRight[0]-fLowLeft[0])/numdomains;
        fMHMSpacing[1] = (fTopRight[1]-fLowLeft[1])/numdomains;
        fMHMSpacingInt[0] = fMHMSpacing[0]/fTol;
        fMHMSpacingInt[1] = fMHMSpacing[1]/fTol;
        fMHMSpacing[0] = fMHMSpacingInt[0]*fTol;
        fMHMSpacing[1] = fMHMSpacingInt[0]*fTol;
    }
    
    void SetTol(REAL tol)
    {
        fTol = tol;
        SetMHMSpacing(fNumMHMDomains);
    }
    /// will project the endnodes on the MHM mesh if they fall within the tolerance
    int64_t InsertNode(TPZGeoNode &gnode);
    
    /// split the fractures if they intersect
    void ComputeFractureIntersections();
    
    /// split the fractures if they intersect
    void ComputeFractureOverlap();
    
    /// add the cornernodes of the MHM mesh
    void AddMHMNodes();
    
    /// return the index of the MHM domain of a fracture
    int MHMDomain(TPZFracture &frac);
    
    /// return an integer pair corresponding to a coordinate
    uint64_t GetLoc(TPZGeoNode &node)
    {
        TPZManVector<REAL, 3> loc(3);
        node.GetCoordinates(loc);
        auto respair = NodeKey(loc);
        
        uint32_t a = respair.first;
        uint32_t b = respair.second;
        uint64_t res = ((uint64_t)a << 32) | b;
        return res;
    }
    
    std::pair<uint32_t,uint32_t> NodeKey(const TPZVec<REAL> &pt)
    {
        TPZManVector<REAL,3> ptcop(pt);
        for(int i=0;i<3;i++) {
            ptcop[i]-=fLowLeft[i];
            ptcop[i] /= fTol;
        }
        std::pair<uint32_t, uint32_t> res((uint32_t)(ptcop[0]+0.5),(uint32_t)(ptcop[1]+0.5));
        return res;
    }
    
    std::pair<uint32_t,uint32_t> NodeKey(int64_t index)
    {
        TPZManVector<REAL,3> co(3,0.);
        fNodeVec[index].GetCoordinates(co);
        return NodeKey(co);
    }
    
    TPZVec<REAL> ConvertNodekey(uint64_t nodekey)
    {
        uint32_t a,b;
        b = nodekey % (1L<<32);
        a = nodekey / (1L<<32);
        TPZManVector<REAL,3> res(3,0.);
        res[0] = fLowLeft[0] + a*fTol;
        res[1] = fLowLeft[1] + b*fTol;
     
        return res;
    }
    
    /// verify the data structure consistency
    void CheckPointMapConsistency();
    
    /// split the fractures
    // if jfrac == -1, split only the ifrac
    void SplitFractures(int64_t ifrac, double xy, int64_t jfrac, double ab);
    
    /// Split fractures by the MHM grid
    void SplitFracturesByMHM();
    
    /// CleanUp FractureNetwork
    void CleanupFractureNetwork();
    
    /// Transfer to the geometric mesh
    void ToGeoMesh();
    
    /// Transfer from the geometric mesh
    void FromGeoMesh();
    
    /// Merge lines which are parallel
    void MergeParallelLines();
    
    /// Delete very short fractures
    void DeleteVeryShortFractures(REAL length);
    
    /// Compute the mesh size at the nodes
    void ComputeMeshSizeAtNodes();
    
};


#endif
