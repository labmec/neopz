#ifndef TPZGEOMESHTOOLS_H
#define TPZGEOMESHTOOLS_H

#include "pzgmesh.h"
#include "MMeshType.h"
/**
 * This namespace concentrates general purpose auxiliary methods for geometric meshes
 */
namespace TPZGeoMeshTools{
    /**
     * This method divides any pyramids (without subelements) in the geometric mesh into two tetrahedra.
     * This method is specially useful in the context of defining composite shape functions in the pyramid
     * based on two tetrahedra and constraints over their sides over the quadrilateral face.
     * @param gmesh
     */
    void DividePyramidsIntoTetra(TPZGeoMesh *gmesh);

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
    TPZGeoMesh * CreateGeoMeshOnGrid(int dim, const TPZVec<REAL>& minX, const TPZVec<REAL>& maxX,
            const TPZVec<int> &matids, const TPZVec<int> nDivs, MMeshType meshType, bool createBoundEls);


    /*! Create one-dimensional geometric mesh.
      \param minX x-coordinate of the leftmost point
      \param maxX x-coordinate of the rightmost point
      \param nEls number of elements to create
      \param matids vector containing the identifier for the materials in the mesh (volume and boundary elements)
      \param createBoundEls whether to create boundary elements
     */
    TPZGeoMesh *
    CreateGeoMesh1D(const REAL minX, const REAL maxX, const int nEls,
                    const TPZVec<int> &matids, bool createBoundEls);

    /*! Create a geometric mesh with a single element corresponding to the reference element.
      It calls TPZGeoMeshTools::CreateGeoMeshSingleElT.
    */
    TPZGeoMesh *
    CreateGeoMeshSingleEl(const MMeshType meshType, const int matid, const bool createBoundEls, const int matidbc = -1);
    /*! Create a geometric mesh with a single element corresponding to the reference element. Optionally creates boundary elements as well.
    */
    template<class TGEO>
    TPZGeoMesh *
    CreateGeoMeshSingleElT(const int matid, const bool createBoundEls, const int matidbc = -1);
}

#endif