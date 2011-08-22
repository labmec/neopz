/**
 * @file
 * @brief Contains the TPZGenGrid class which implements the generation of a multilayered geometric grid.
 */
//
// Author: MISAEL LUIS SANTANA MANDUJANO.
//
// File:   TGenGrid.h
//
// Class:  TGenGrid
//
// Obs.:   Gera uma malha retangular:
//
// Versao: 06 / 1996.
//

#ifndef _TPZGENGRIDHH_
#define _TPZGENGRIDHH_

class TPZCompMesh;
class TPZGeoMesh;
#include "pzvec.h"

#include <stdio.h>
#include <iostream>
#include "pzreal.h"
#include "tpzautopointer.h"

#include <fstream>

/** 
 * @ingroup pre
 * @brief Implements the generation of a multilayered geometric grid. \ref pre "Getting Data"
 */
/** This class uses DEPRECATED objects, but can be easily updated
 */
class TPZGenGrid{
    
public:
    
    /**
	 @brief Constructor
	 @param nx 
     @param x0 lower left coordinate
     @param x1 upper right coordinate
     @param numl number of layers
     @param rot rotation applied to the grid
     */
    TPZGenGrid(TPZVec<int> &nx, TPZVec<REAL> &x0, TPZVec<REAL> &x1, int numl = 1, REAL rot = 0.5);
    
    virtual ~TPZGenGrid();
    
    /** @brief Add nodes and elements to the object mesh */
    virtual short Read (TPZGeoMesh & malha);
    
    /** @brief Compute the geometric progression such that the first elements have this size */
    static REAL GeometricProgression(REAL minsize, REAL size, int numdiv);
	
    /** @brief Compute the geometric progression such that the first elements have this size */
    void ComputeGeometricProgression(TPZVec<REAL> &minsizes, TPZVec<REAL> &progression);
    
    /** @brief Sets the geometric progression of the mesh to be generated */
    void SetGeometricProgression(TPZVec<REAL> &progression);
    
    /** @brief Generate boundary geometric elements associated with the side
     * @note Note : range of the side parameter 0 <= side < 4
     */
    virtual void SetBC(TPZGeoMesh *gr, int side, int bc);
    
    /** @brief Generate boundary geometric elements between node start and end going counter clockwise */
    virtual void SetBC(TPZGeoMesh *g, TPZVec<REAL> &start, TPZVec<REAL> &end, int bc);
    
    /** @brief Generate boundary geometric elements associated with the side
     * @note Note : range of the side parameter 0 <= side < 4
     */
    void SetBC(TPZAutoPointer<TPZGeoMesh> gr, int side, int bc)
    {
        SetBC(gr.operator->(), side, bc);
    }
    
    /** @brief Generate boundary geometric elements between node start and end going counter clockwise */
    void SetBC(TPZAutoPointer<TPZGeoMesh> gr, TPZVec<REAL> &start, TPZVec<REAL> &end, int bc)
    {
        SetBC(gr.operator->(), start, end, bc);
    }
    
    /** @brief Generate a boundary geometric element at the indicated node */
    void SetPointBC(TPZGeoMesh *gr, TPZVec<REAL> &x, int bc);
	
    /** @brief Generate a boundary geometric element at the indicated node */
    void SetPointBC(TPZAutoPointer<TPZGeoMesh> gr, TPZVec<REAL> &x, int bc)
    {
        SetPointBC(gr.operator->(), x, bc);
    }
	
    /** @brief Prints the data structure of the class */
    virtual void Print( char *name = NULL, std::ostream &out = std::cout );
    
    /** @brief Set the element type
	 *
	 * type = 0 -> quadrilateral
	 * type = 1 -> triangle
	 * type = 2 -> quadratic quadrilaterals
     */
    virtual void SetElementType(int type);
    
    /** @brief Returns the element id for the element addressed by the parameters */
    int ElemId(int iel,int jel, int layer);
    
    /** @brief Computes the distance between two points */
	static REAL Distance(TPZVec<REAL> &x1,TPZVec<REAL> &x2);
	
protected:
    
    virtual void Coord(int i, TPZVec<REAL> &coord);
    
    virtual int GlobalI(int ix, int iy, int layer);
    
	void ElementConnectivity(int iel, TPZVec<int> &nodes);
    
    virtual void GenerateNodes(TPZGeoMesh &grid);
    
    virtual void GenerateElements(TPZGeoMesh &grid);
    
    /** @brief Number of elements in both directions */
	TPZVec<int> fNx;
    /** @brief Coordinate of the lower left point */
	TPZVec<REAL> fX0;
    /** @brief coordinate of the upper right point */
	TPZVec<REAL> fX1;
    /** @brief Size of the lower left element */
	TPZVec<REAL> fDelx;
    /** @brief Geometric progression coeficients in the x and y direction */
	TPZVec<REAL> fGeometricProgression;
    /** @brief Number of nodes of the mesh */
	int fNumNodes;
    /** 
	 * @brief Variable which indicates the type of element that should be generated
     * type = 0 -> quadrilateral
     * type = 1 -> triangle
     * type = 2 -> quadratic quadrilaterals
     */
	int fElementType;
    
    /**
     @brief Number of meshes which will be generated hinging along an axis
     */
    int fNumLayers;
    /** 
     @brief Rotation angle between the layers
     */
    REAL fRotAngle;
	
};

#endif // _TGENGRIDHH_
