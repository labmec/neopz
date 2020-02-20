//
// Created by Francisco Teixeira Orlandini on 19/12/19.
//

#ifndef PZ_TPZGENGRID3D_H
#define PZ_TPZGENGRID3D_H
#include "pzreal.h"
#include "pzeltype.h"
#include "MMeshType.h"
#include "pzmanvector.h"
class TPZGeoMesh;


/**
 * This class generates a uniform three-dimensional mesh with tetrahedral, hexahedral or prismatic
 * elements. It could be extended with little effort for generating pyramidal elements as well.
 */
class TPZGenGrid3D {
public:
    TPZGenGrid3D() = delete;

    /**
     * Initializes the structure for creating a mesh that forms a partition of a cuboid.
     * minX and maxX are spatial coordinates of the endpoints of one spatial diagonal of the cuboid.
     * The cuboid faces are aligned with the x,y and z planes.
     * The mesh is structured and its divisions are given by nelDiv.
     * @param minX The smallest (or most negative) values of x,y and z for a point inside the cuboid.
     * @param maxX The biggest (or least negative) values of x,y and z for a point inside the cuboid.
     * @param nelDiv Number of elements in each direction (x,y and z)
     * @param elType How to divide each cuboid of the mesh.
     */
    TPZGenGrid3D(const TPZVec<REAL> &minX, const TPZVec<REAL> &maxX, const TPZVec<int> &nelDiv, const MMeshType elType);

    /**
     * This method builds the three-dimensional elements of the desired type and sets their material ids as matIdDomain.
     * The generated elements are of type TPZGeoElRefPattern<elType> thus can be refined if desired.
     * @param matIdDomain Material identifier of the volumetric elements
     */
    TPZGeoMesh *BuildVolumetricElements(const int matIdDomain);

    /**
     * This method creates the boundary elements and associates them with their respective material ids.
     * Their order follows the enumeration of the sides of a hexahedral element.
     * @param matIdZmin Material identifier for the boundary elements on the min z coordinate face
     * @param matIdXmin Material identifier for the boundary elements on the min x coordinate face
     * @param matIdYmin Material identifier for the boundary elements on the min y coordinate face
     * @param matIdXmax Material identifier for the boundary elements on the max x coordinate face
     * @param matIdYmax Material identifier for the boundary elements on the max y coordinate face
     * @param matIdZmax Material identifier for the boundary elements on the max z coordinate face
     * @return
     */
    TPZGeoMesh *BuildBoundaryElements(const int matIdZmin, const int matIdXmin,
                                      const int matIdYmin, const int matIdXmax,
                                      const int matIdYmax, const int matIdZmax);
private:
    TPZGeoMesh *fGmesh{nullptr};
    const MMeshType fMeshType;
    const TPZManVector<int,3>fNelDiv{0};
    const TPZManVector<REAL,3> fMinX{0};
    const TPZManVector<REAL,3> fMaxX{0};
    static const int fDim{3};
};


#endif //PZ_TPZGENGRID3D_H
