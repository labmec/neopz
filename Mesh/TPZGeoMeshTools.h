#ifndef TPZGEOMESHTOOLS_H
#define TPZGEOMESHTOOLS_H

#include "pzgmesh.h"
#include "MMeshType.h"
/**
 * This class provides general purpose auxiliary methods for geometric meshes
 */
class TPZGeoMeshTools{
public:

    /**
     * This method divides any pyramids (without subelements) in the geometric mesh into two tetrahedra.
     * This method is specially useful in the context of defining composite shape functions in the pyramid
     * based on two tetrahedra and constraints over their sides over the quadrilateral face.
     * @param gmesh
     */
    static void DividePyramidsIntoTetra(TPZGeoMesh *gmesh);

    /**
     * This method creates a TPZGeoMesh of dim 3(2) based on a cuboid(rectangle) aligned with the x,y and z axis (x and y).
     * @param dim Mesh dimension
     * @param minX The smallest (or most negative) values of x,y and z for a point inside the cuboid (rectangle).
     * @param maxX The biggest (or least negative) values of x,y and z for a point inside the rectangle.
     * @param matids Material identifiers. If creating boundaries, it must have size dim*2 + 1.
     * Their order follows the side numbering of the hexahedral element (quadrilateral element)
     * @param nDivs Number of elements in each direction
     * @param meshType What type of elements to create
     * @param createBoundEls Whether to create the boundary elements
     * @return Pointer to the TPZGeoMesh that has been created
     */
    static TPZGeoMesh * CreateGeoMeshOnGrid(int dim, const TPZVec<REAL>& minX, const TPZVec<REAL>& maxX,
            const TPZVec<int> &matids, const TPZVec<int> nDivs, MMeshType meshType, bool createBoundEls);
};

#endif