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

public:
TPZHierarquicalGrid();
TPZHierarquicalGrid(TPZGeoMesh *Geomesh);
TPZHierarquicalGrid(const TPZHierarquicalGrid& other);
virtual ~TPZHierarquicalGrid();
virtual TPZHierarquicalGrid& operator=(const TPZHierarquicalGrid& other);
virtual bool operator==(const TPZHierarquicalGrid& other) const;

    /**
     * @brief Thickness of the mesh (+  or -)
     */
    REAL fThickness;
    
    /**
     * @brief Name of the fine mesh to be extended
     */
    std::string fFileName;
    
    /**
     * @brief A geometric mesh generated from other sources
     */
    TPZAutoPointer<TPZGeoMesh> fBase;
    
    /**
     * @brief Vector of n bases to be connected to the original base
     */
    TPZVec<TPZAutoPointer<TPZGeoMesh> > fSubBases;
    
    /** @brief Pointer to parametric function of t parameter */
    TPZAutoPointer<TPZFunction<STATE> > fParametricFunction;    
    
    /**
     * @brief Prints the generated mesh
     */
    void PrintGeneratedMesh(std::ostream &out = std::cout);    

    // Set Get Methods
    
    void SetParametricFunction(TPZAutoPointer<TPZFunction<STATE> > fp)
    {
        fParametricFunction = fp;
    }
    
    TPZGeoMesh * ComputeExtrusion(STATE t, STATE dt, int n);
    
    void CreateGeometricElement(int n, int iel, int eldim, int elmatid, int &elid, TPZGeoMesh * ngmesh);
    
};

#endif // TPZHIERARQUICALGRID_H
