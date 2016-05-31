//
//  TPZAcademicGeoMesh.hpp
//  PZ
//
//  Created by Philippe Devloo on 5/29/16.
//
//

#ifndef TPZAcademicGeoMesh_hpp
#define TPZAcademicGeoMesh_hpp

#include <stdio.h>
#include "pzmanvector.h"
#include "pzreal.h"
#include "pzgmesh.h"



class TPZAcademicGeoMesh
{
public:
    static TPZManVector<REAL,3> fX0, fEps;
    
    TPZGeoMesh fDeformed;
    
    enum MMeshType {ENone, EHexa, EPrism, ETetrahedra, EPyramid };

protected:
    
    /// type of mesh that should be generated
    MMeshType fMeshType;
    
    /// indices of the boundary conditions
    TPZManVector<int,6> fBCNumbers;
    
    /// whether the mesh should be deformed or not
    bool fShouldDeform;
    
    /// number of elements in any direction
    long fNumberElements;
    
    /// material index for volume elements
    int fMaterialId;
    
    /// put the geometric nodes in the geometric mesh
    void GenerateNodes(TPZGeoMesh *gmesh);
    
    /// Deform the geometric mesh according to the coordinates of fDeformed
    void DeformGMesh(TPZGeoMesh &gmesh);
    
    /** Constructing geometrical mesh depends on type of element wished. */
    TPZGeoMesh *HexahedralMesh();
    TPZGeoMesh *PyramidalAndTetrahedralMesh();
    TPZGeoMesh *TetrahedralMesh();
    TPZGeoMesh *TetrahedralMeshUsingRefinement();
    
    int AddBoundaryElements(TPZGeoMesh *gmesh);
    
    void CheckConsistency(TPZGeoMesh *gmesh);

    
public:
    
    /// sets up parameters to create de hexahedral regular mesh
    TPZAcademicGeoMesh();
    
    /// sets up parameters to create de hexahedral regular mesh
    TPZAcademicGeoMesh(int numelements, MMeshType meshtype, bool deform = false);
    
    /// set the meshtype acording to the enumerate values
    void SetMeshType(MMeshType meshtype)
    {
        if (meshtype == ENone) {
            DebugStop();
        }
        fMeshType = meshtype;
    }
    
    /// set the boundary numbers on the sides of the cube
    void SetBoundaryIndices(TPZVec<int> &BCNumbers);
    
    /// toggle between regular or deformed mesh
    void ShouldDeform(bool deform);
    
    /// set the number of elements in any direction
    void SetNumberElements(long numelements)
    {
#ifdef PZDEBUG
        if(numelements < 0)
        {
            DebugStop();
        }
#endif
        fNumberElements = numelements;
    }
    /// set the material id to be used
    void SetMaterialId(int id)
    {
#ifdef PZDEBUG
        if(id == 0)
        {
            DebugStop();
        }
#endif
        fMaterialId = id;
    }
    
    /// create a geometric mesh acording to the parameters of the object
    TPZGeoMesh *CreateGeoMesh();
    
    
};

#endif /* TPZAcademicGeoMesh_hpp */
