/**
 * @file
 * @brief Contains the TPZGenGrid class which implements the generation of a multilayered geometric grid (two-dimensional).
 */
//
// Author: MISAEL LUIS SANTANA MANDUJANO.
//
// File:   TGenGrid.h
//
// Class:  TGenGrid
//
// Obs.:   Gera uma malha bi-dimensional:
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
 * @brief Implements the generation of a multilayered bi-dimensional geometric grid. \ref pre "Getting Data"
 */
/** This class uses DEPRECATED objects, but can be easily updated */
class TPZGenGrid {

public:

    /**
	 * @brief Constructor of the rectangular domain
	 * @param nx numbers of partition intervals for x (nx[0]) and for y (nx[1])
     * @param x0 lower left coordinate
     * @param x1 upper right coordinate
     * @param numl number of layers
     * @param rot rotation applied to the grid for the next layer
	 * @note All the layers has a common interval \f$ [(x0[0],0,0);(x1[0],0,0)] \f$
     */
    TPZGenGrid(TPZVec<int> &nx, TPZVec<REAL> &x0, TPZVec<REAL> &x1, int numl = 1, REAL rot = 0.5);

	/** @brief Default destructor */
    virtual ~TPZGenGrid();

    /**
	 * @brief Add nodes and elements to the object mesh
	 * @param mesh Object mesh for which will be created the nodes and elements (depends on fTypeElement)
	 */
    virtual short Read(TPZGeoMesh &mesh);
    /**
	 * @brief Add nodes and elements to the object mesh
	 * @param mesh Object mesh for which will be created the nodes and elements (depends on fTypeElement)
	 */
    short Read(TPZAutoPointer<TPZGeoMesh> mesh);
    
    /**
	 * @brief Compute the geometric progression such that the first elements have this size
	 */
    static REAL GeometricProgression(REAL minsize, REAL size, int numdiv);
	
    /** @brief Compute the geometric progression such that the first elements have this size */
    void ComputeGeometricProgression(TPZVec<REAL> &minsizes, TPZVec<REAL> &progression);
    
    /** @brief Sets the geometric progression of the mesh to be generated */
    void SetGeometricProgression(TPZVec<REAL> &progression);
    
    /**
	 * @brief Generate boundary geometric elements associated with the side of the rectangular domain
	 * @param gr object mesh
	 * @param side boundary side of the rectangular domain which will be set the boundary condition
	 * @param bc boundary condition id (material), generally negative
     * @note Note : range of the side parameter 0 <= side < 4
     */
    virtual void SetBC(TPZGeoMesh *gr, int side, int bc);
    
    /**
	 * @brief Generate boundary geometric elements between node start and end going counter clockwise
	 * @param gr Object mesh which will be set the elements with the specified boundary condition
	 * @param start First node (boundary coordinates) which will be identified with boundary id
	 * @param end Last node (boundary coordinates) which will identified with boundary id
	 * @param bc Boundary condition id (material), generally negative
	 * @note From start to end going counter clockwise
	 */
    virtual void SetBC(TPZGeoMesh *gr, TPZVec<REAL> &start, TPZVec<REAL> &end, int bc);
    
    /**
	 * @brief Generate boundary geometric elements associated with the side of the rectangular domain
	 * @param gr object mesh
	 * @param side boundary side of the rectangular domain which will be set the boundary condition
	 * @param bc boundary condition id (material), generally negative
     * @note Note : range of the side parameter 0 <= side < 4
     */
    void SetBC(TPZAutoPointer<TPZGeoMesh> gr, int side, int bc)
    {
        SetBC(gr.operator->(), side, bc);
    }
    
    /**
	 * @brief Generate boundary geometric elements between node start and end going counter clockwise
	 * @param gr Object mesh which will be set the elements with the specified boundary condition
	 * @param start First node (boundary coordinates) which will be identified with boundary id
	 * @param end Last node (boundary coordinates) which will identified with boundary id
	 * @param bc Boundary condition id (material), generally negative
	 */
    void SetBC(TPZAutoPointer<TPZGeoMesh> gr, TPZVec<REAL> &start, TPZVec<REAL> &end, int bc)
    {
        SetBC(gr.operator->(), start, end, bc);
    }
    
    /**
	 * @brief Generate a boundary geometric element at the indicated node
	 * @param gr Object mesh which will be set the elements with the specified boundary condition
	 * @param x Node (boundary coordinates) which will identified with boundary id
	 * @param bc Boundary condition id (material), generally negative
	 */
    void SetPointBC(TPZGeoMesh *gr, TPZVec<REAL> &x, int bc);
	
    /**
	 * @brief Generate a boundary geometric element at the indicated node
	 * @param gr Object mesh which will be set the elements with the specified boundary condition
	 * @param x Node (boundary coordinates) which will identified with boundary id
	 * @param bc Boundary condition id (material), generally negative
	 */
    void SetPointBC(TPZAutoPointer<TPZGeoMesh> gr, TPZVec<REAL> &x, int bc)
    {
        SetPointBC(gr.operator->(), x, bc);
    }
	
    /**
	 * @brief Prints the data structure of the class
	 * @param name Title of the structure data will be printed
	 * @param out Output stream (ofstream of the output file data)
	 */
    virtual void Print( char *name = NULL, std::ostream &out = std::cout );
    
    /**
	 * @brief Set the element type
	 * @param type Bi-dimensional element type
	 * @note Values: \f$ type = 0 \f$ (quadrilateral); \f$ type = 1 \f$ (triangle); \f$ type = 2 \f$ (quadratic quadrilaterals)
     */
    virtual void SetElementType(int type);
    
    /**
	 * @brief Returns the element id for the element addressed by the parameters
	 * @param iel ith division for x
	 * @param jel jth division for y
	 * @param layer number of layer for which the element belongs
	 */
    int ElemId(int iel,int jel, int layer);
    
    /**
	 * @brief Computes the euclidean distance between two points, or the measure of the interval between two points
	 * @param x1 start point
	 * @param x2 last point
	 */
	static REAL Distance(TPZVec<REAL> &x1,TPZVec<REAL> &x2);
	
	/**
	 * @brief Merges two geometrical mesh created for TPZGenGrid as separated
	 * @param grid Mesh over which will be increment the nodes and elements no duplicated of the second mesh
	 * @param grid2 Mesh from get nodes and elements and put into grid if it is not duplicated
	 */
	bool ReadAndMergeGeoMesh(TPZAutoPointer<TPZGeoMesh> grid,TPZAutoPointer<TPZGeoMesh> grid2);
	
protected:
    /**
	 * @brief Computes the coordinates of the ith geometric node generated
	 * @param i generated geometric node id
	 * @param coord coordinates of the ith node
	 */
    virtual void Coord(int i, TPZVec<REAL> &coord);
    
    /**
	 * @brief Returns the geometric node id for the element addressed by the parameters
	 * @param ix ith division for x
	 * @param iy ith division for y
	 * @param layer number of the layer for which the node belongs
	 */
    virtual int GlobalI(int ix, int iy, int layer);
    
	void ElementConnectivity(int iel, TPZVec<int> &nodes);
    
	/**
	 * @brief Creates the geometric nodes, it depends on fElementType, layer and fRotAngle
	 * @param grid Object mesh in which will be stored the created nodes
	 */
    virtual bool GenerateNodes(TPZGeoMesh &grid);
    
	/**
	 * @brief Creates the geometric element: triangles or quadrilaterals
	 * @param grid Object mesh in which will be stored the created elements
	 */
    virtual bool GenerateElements(TPZGeoMesh &grid);
    
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
     * \li type = 0 -> quadrilateral
	 * \li type = 1 -> triangle
     * \li type = 2 -> quadratic quadrilaterals
     */
	int fElementType;
    
    /** @brief Number of meshes which will be generated hinging along an axis */
    int fNumLayers;
    /** @brief Rotation angle between the layers */
    REAL fRotAngle;
	
};

#endif // _TGENGRIDHH_
