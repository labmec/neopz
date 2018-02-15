//
//  TRMGmshReader.h
//  PZ
//
//  Created by Omar on 1/15/17.
//
//

#ifndef TRMGmshReader_h
#define TRMGmshReader_h

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include "pzgmesh.h"

class TPZGeoMesh;


struct MaterialDataS {

    TPZStack<int> fMatID;
    TPZStack<std::pair<int ,std::string> >  fMaterial;
    
    MaterialDataS() : fMatID(), fMaterial(){
        
    }
    
    MaterialDataS(int num) : fMatID(), fMaterial(){
        
    }
    
    MaterialDataS(const MaterialDataS &copy) : fMatID(copy.fMatID),
    fMaterial(copy.fMaterial) {
    }
    
    MaterialDataS &operator=(const MaterialDataS &copy){
        fMatID = copy.fMatID;
        fMaterial = copy.fMaterial;
        return *this;
    }
    
};


/**
 * @brief Implement the interface between TPZGeoMesh and the files produced by Gmsh (version 2.2 ) in msh format.
 * @since January 16, 2017
 */

/** What is Gmsh ? Take a look on http://gmsh.info/
 * Gmsh is a free 3D finite element grid generator with a build-in CAD engine and post-processor. Its design goal is to provide a fast, light and user-friendly
 * meshing tool with parametric input and advanced visualization capabilities. Gmsh is built around four modules: geometry, mesh, solver and post-processing.
 * The specification of any input to these modules is done either interactively using the graphical user interface or in ASCII text files using Gmsh's own
 * scripting language.
 */
class TRMGmshReader{
    
public:

    /** @brief default destructor */
    TRMGmshReader();
    
    /** @brief default destructor */
    ~TRMGmshReader();
    
    /** @brief Convert Gmsh msh files in a TPZGeoMesh object */
    TPZGeoMesh * GeometricGmshMesh(std::string file_name);
    
    /** @brief Number of Materials */
    /** Number of volumetric materials */
    int fVolNumber;
    
    /** @brief Number of Boundary Conditions */
    /** Number of Boundary Conditions */
    int fBCNumber;
    
    /** @brief Mesh Dimension */
    /** Mesh Dimension */
    int fProblemDimension;
    
    /** @brief Mesh Dimension */
    /** Mesh Dimension */
    REAL fDimensionlessL;
    
    /** @brief MaterialVec */
    /** Structure of both, physical entities dimension and names */
    MaterialDataS fMaterialDataVec;
    
    
    /** @brief Characteristic domain dimension for dimensionless geometry. */
    /** Set Max dimension for geometric domain default = 1.0. */
    void SetfDimensionlessL(REAL dimensionlessL);
    
    /** @brief Insert elements following msh file format */
    bool InsertElement(TPZGeoMesh * gmesh, std::ifstream & line);
    
private:
    
    
};

#endif /* TRMGmshReader_h */
