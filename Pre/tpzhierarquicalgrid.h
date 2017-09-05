/*
 <one line to give the program's name and a brief idea of what it does.>
 Copyright (C) 2014  <copyright holder> <email>
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#ifndef TPZHIERARQUICALGRID_H
#define TPZHIERARQUICALGRID_H

class TPZGeoMesh;
#include "pzstack.h"
#include <fstream>
#include <iostream>
#include "pzmatrix.h"
#include "tpzautopointer.h"
#include "pzfunction.h"

class TPZHierarquicalGrid
{

    /**
     * @brief Thickness of the mesh (+  or -)
     */
    REAL fThickness;
    
    /**
     * @brief Non affine 2D and 3D extrusion
     */
    bool fNonAffineQ;
    
    /**
     * @brief 2d extrusion is based on quadrilaterals
     */
    bool fIsQuad;
    
    /**
     * @brief 3d extrusion is based on prisms
     */
    bool fIsPrism;
    
    /**
     * @brief 3d extrusion is based on tetrahedrons
     */
    bool fIsTetrahedron;
    
    /**
     * @brief Extrusion front material id
     */
    int ffrontMatID;
    
    /**
     * @brief Extrusion back material id
     */
    int fbackMatID;
    
    /**
     * @brief Name of the fine mesh to be extended
     */
    std::string fFileName;
    
    /**
     * @brief A geometric mesh being computed
     */
    TPZGeoMesh * fComputedGeomesh;
    
    /**
     * @brief A geometric mesh generated from other sources
     */
    TPZAutoPointer<TPZGeoMesh> fBase;
    
    /**
     * @brief Vector of n bases to be connected to the original base
     */
    TPZVec<TPZAutoPointer<TPZGeoMesh> > fSubBases;
    
    /** @brief Pointer to parametric function of t parameter */
    TPZAutoPointer<TPZFunction<REAL> > fParametricFunction;
    

public:
    
    /**
     *  Default constructor
     */
    TPZHierarquicalGrid();
    
    /**
     *  TPZGeoMesh based constructor
     */
    TPZHierarquicalGrid(TPZGeoMesh *Geomesh);
    
    /**
     *  Copy constructor
     */
    TPZHierarquicalGrid(const TPZHierarquicalGrid& other);
    
    /**
     *  Destructor
     */
    ~TPZHierarquicalGrid();
    
    /**
     *  Assignment operator
     *
     *  @param other TPZHierarquicalGrid object being copied
     *
     *  @return The new copy
     */
    TPZHierarquicalGrid& operator=(const TPZHierarquicalGrid& other);
    
    /**
     *  Comparison operator
     *
     *  @param other right hand side
     *
     *  @return comparison result (true/false)
     */
    bool operator==(const TPZHierarquicalGrid& other) const;
    
    /**
     * @brief Prints the generated mesh
     */
    void PrintGeneratedMesh(std::ostream &out = std::cout);
    
    /**
     *  Set the geometric mesh being extruded
     *
     *  @param gmesh TPZGeoMesh representation of the geometry
     */
    void SetGeometricMesh(TPZGeoMesh * gmesh){
        fBase = gmesh;
    }
    
    /**
     *  Paramatric function used during the extrusion
     *
     *  @param fp autopointer  of the parametric function
     */
    void SetParametricFunction(TPZAutoPointer<TPZFunction<REAL> > fp)
    {
        fParametricFunction = fp;
    }
    
    /**
     *  Compute the extrusion based on a lower dimension mesh
     *
     *  @param t  one-dimensional parameteric representation of the R3 space
     *  @param dt element size over the parametric space
     *  @param n  number of elements in the paremetric direction
     *
     *  @return high dimension geometric mesh
     */
    TPZGeoMesh * ComputeExtrusion(REAL t, REAL dt, int n);
    
    void SetFrontBackMatId(int front, int back) {ffrontMatID = front; fbackMatID = back;}
    
    void SetTriangleExtrusion() {fIsQuad = false;}
    
    void SetTetrahedonExtrusion() {fIsTetrahedron = true;}
    
    void SetPrismExtrusion() {fIsPrism = true;}
    
    void SetNonAffineExtrusion() { fNonAffineQ = true;}
    
    void SetGridFileName(std::string &FileName) {fFileName = FileName;}
    
    void CreateGeometricElement(int n, int iel, int eldim, int elmatid, int &elid);
    
};

#endif // TPZHIERARQUICALGRID_H
