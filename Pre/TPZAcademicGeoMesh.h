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
#include <iostream>


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
    int64_t fNumberElements;
    
    /// material index for volume elements
    int fMaterialId;
    
    /// put the geometric nodes in the geometric mesh
    void GenerateNodes(TPZGeoMesh *gmesh);
    
    /// Deform the geometric mesh according to the coordinates of fDeformed
    void DeformGMesh(TPZGeoMesh &gmesh);
    
    /** Constructing geometrical mesh depends on type of element wished. */
    TPZGeoMesh *HexahedralMesh();
    TPZGeoMesh *TetrahedralMeshUsingRefinement();
    
    int AddBoundaryElements(TPZGeoMesh *gmesh);
  
    int AddBoundaryElementsByCoord(TPZGeoMesh *gmesh);
      
    void CheckConsistency(TPZGeoMesh *gmesh);

    
public:
    
    /// sets up parameters to create de hexahedral regular mesh
    TPZAcademicGeoMesh();
    
    /// sets up parameters to create de hexahedral regular mesh
    TPZAcademicGeoMesh(int numelements, MMeshType meshtype, bool deform = false) : fMeshType(meshtype), fBCNumbers(6,-1), fShouldDeform(deform),
    fNumberElements(numelements), fMaterialId(1)
    {
        
    }
    
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
    void SetNumberElements(int64_t numelements)
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
    TPZGeoMesh *CreateGeoMesh()
    {
        switch (fMeshType) {
            case EHexa:
                return HexahedralMesh();
                break;
            case ETetrahedra:
                return TetrahedralMesh();
                break;
            case EPyramid:
                return PyramidalAndTetrahedralMesh();
                break;
            case EPrism:
                std::cout << "Not implemented\n";
                
            default:
                DebugStop();
                break;
        }
        return 0;
    }
  
    TPZGeoMesh *PyramidalAndTetrahedralMesh();
    TPZGeoMesh *RedBlackPyramidalAndHexagonalMesh();
    TPZGeoMesh *TetrahedralMesh();
  
    /**
    *  Sets the bc id vector for the hexaedral mesh
    *
    *  @param BCNumbers id of bcs with order and location on global mesh
    *  0: x = 0
    *  1: y = 0
    *  2: z = 0
    *  3: x = 1
    *  4: y = 1
    *  5: z = 1
    */
  void SetBCIDVector(TPZVec<int> &BCNumbers){
    if(fBCNumbers.size() != BCNumbers.size()){
      DebugStop();
    }
    for (int i = 0; i < BCNumbers.size(); i++) {
      fBCNumbers[i] = BCNumbers[i];
    }
  }
  
    
};

#endif /* TPZAcademicGeoMesh_hpp */
