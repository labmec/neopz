#ifndef TPZGEOMESHTOOLS_H
#define TPZGEOMESHTOOLS_H

#include "pzgmesh.h"
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
};

#endif