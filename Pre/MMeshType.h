#ifndef MMESHTYPE_H
#define MMESHTYPE_H

/**
 * @brief This class defined mesh types for usage of mesh generating classes TPZGenGrid2D and TPZGenGrid3D.
 */
enum class MMeshType
{
    EQuadrilateral,/*For 2D meshes consisting of quadrilaterals*/
    ETriangular,/*For 2D meshes where each quadrilateral is split into two triangles*/
    EHexahedral,/*For 3D meshes consisting of hexahedra*/
    ETetrahedral,/*For 3D meshes where each hexahedron is divided into five tetrahedra*/
    EPyramidal,/*For 3D meshes where each hexahedron is divided into six pyramids*/
    EPrismatic,/*For 3D meshes where each hexahedron is divided into two prisms*/
    EHexaPyrMixed,/*For 3D meshes where alternating cubes are divided into six pyramids*/
    ENoType
};

/**
 * Returns the name associated with a certain mesh type
 * @param meshType mesh type
 * @return name
 */
std::string MMeshType_Name(const MMeshType meshType);

/**
 * Returns the dimension of the mesh associated with a given mesh type
 * @param meshType mesh type
 * @return dimension
 */
int MMeshType_Dimension(const MMeshType meshType);

/**
 * Stream operator for getting the mesh type name
 * @param out
 * @param meshType
 * @return
 */
std::ostream &operator<<(std::ostream &out, const MMeshType meshType);

#endif