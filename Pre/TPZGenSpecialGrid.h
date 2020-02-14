/**
 * @file
 * @brief Declaration of the generator of special grids class.
 * @note The idea is to generate a polygonal mesh approximating a three dimensional geometric element from another simple polygonal mesh.
 */

#ifndef GENERATINGSPECIALGRIDSHH
#define GENERATINGSPECIALGRIDSHH

#include "pzreal.h"

#include "pzgmesh.h"

/** 
 * @ingroup pre
 * @brief Implements the generation of a polygonal mesh approximating a geometric element from another simple polygonal mesh. \ref pre "Getting Data"
 */
class TPZGenSpecialGrid {

public:
	
	/** 
	 * @brief Static function to generate a polygonal mesh approximating a sphere from a octahedron mesh (polygonal) 
	 * @param center Center of the sphere
	 * @param radius Radius of the sphere polygonalized
	 * @param nUniformRefs Number of uniform refinements to make.
	 */
	static TPZGeoMesh *GeneratePolygonalSphereFromOctahedron(TPZVec<REAL> &center,REAL radius, int nUniformRefs);

	/** 
	 * @brief Static function to generate a polygonal mesh approximating a sphere from a octahedron mesh (polygonal) 
	 * @param center Center of the sphere
	 * @param radius Radius of the sphere polygonalized
	 * @param tol Tolerance of approximation. It is the distance of the barycenter of any triangular element of the mesh until the sphere.
	 * @note tol is the stop criterion.
	 */
	static TPZGeoMesh *GeneratePolygonalSphereFromOctahedron(TPZVec<REAL> &center,REAL radius, REAL tol);
 
	/** 
	 * @brief Make uniform refinement of the geometric mesh. This method finishs making the new connectivities for the elements (ResetConnectivities() and BuildConnectivity()) 
	 * @param nUniformRefs Number of divisions of all elements wished
	 * @param gmesh Geometric mesh contained the elements
	 * @param dimelfordivision Dimension of the geometric elements will be divided. If dimelfordivision is negative subdivide all elements (for any dimensions)
	 * @param allmaterial If true divides all geometric element, if false check the material match with allmaterial to divide the element
	 * @param matidtodivided To divide only elements with this material id
	 */
	static void UniformRefinement(const int nUniformRefs, TPZGeoMesh *gmesh, const int dimelfordivision,bool allmaterial=true, const int matidtodivided=0);
	
};

#endif