//
// Created by Francisco Teixeira Orlandini on 19/12/19.
//

#ifndef PZ_TPZGENGRID3D_H
#define PZ_TPZGENGRID3D_H
#include "pzreal.h"
#include "pzeltype.h"

class TPZGeoMesh;

enum class MMeshType
{
    ETetrahedral,/*each cube is divided into five tetrahedra*/
    EPyramidal,/*each cube is divided into six pyramids*/
    EPrismatic,/*each cube is divided into two prisms*/
    EHexahedral,
    EHexaPyrMixed,/*alternating cubes are divided into six pyramids*/
    ENoType
};

std::ostream &operator<<(std::ostream &out, const MMeshType meshType);

/**
 * This class generates a uniform three-dimensional mesh with tetrahedral, hexahedral or prismatic
 * elements. It could be extended with little effort for generating pyramidal elements as well.
 */
class TPZGenGrid3D {
public:
    TPZGenGrid3D() = delete;
    TPZGenGrid3D(const REAL minX, const REAL minY, const REAL minZ,const REAL maxX, const REAL maxY, const REAL maxZ,
            const int nelx, const int nely, const int nelz, const MMeshType elType);

    /**
     * This method builds the three-dimensional elements of the desired type and sets their material ids as matIdDomain.
     * The generated elements are of type TPZGeoElRefPattern<elType> thus can be refined if desired.
     * @param matIdDomain Material identifier of the volumetric elements
     */
    TPZGeoMesh *BuildVolumetricElements(const int matIdDomain);

    /**
     * This method creates the boundary elements and associates them with their respective material ids
     * @param matIdXmin Material identifier for the boundary elements on the min x coordinate face
     * @param matIdXmax Material identifier for the boundary elements on the max x coordinate face
     * @param matIdYmin Material identifier for the boundary elements on the min y coordinate face
     * @param matIdYmax Material identifier for the boundary elements on the max y coordinate face
     * @param matIdZmin Material identifier for the boundary elements on the min z coordinate face
     * @param matIdZmax Material identifier for the boundary elements on the max z coordinate face
     * @return
     */
    TPZGeoMesh *BuildBoundaryElements(const int matIdXmin, const int matIdXmax,
                                      const int matIdYmin, const int matIdYmax,
                                      const int matIdZmin, const int matIdZmax);
private:
    TPZGeoMesh *fGmesh{nullptr};
    const MMeshType fMeshType{MMeshType::ENoType};
    const int fNelx{0},fNely{0},fNelz{0};
    const REAL fMinX{0},fMinY{0},fMinZ{0};
    const REAL fMaxX{0},fMaxY{0},fMaxZ{0};
    static const int fDim{3};
};


#endif //PZ_TPZGENGRID3D_H
