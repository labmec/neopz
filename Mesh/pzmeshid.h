/**
 * @file
 * @brief Contains the Ids for several geometric structures.
 */
#ifndef PZMESHIDH
#define PZMESHIDH
//
// C++ Interface: pzmeshid
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@corona>, (C) 2004
//
// Copyright: See COPYING file that comes with this distribution
//
//

/** \addtogroup geometry
 * @{
 */
/** @brief Id of geometric point element */
const int TPZFGEOELEMENTPOINTID = 200;
/** @brief Id of geometric linear element */
const int TPZFGEOELEMENTLINEARID = 201;
/** @brief Id of geometric triangular element */
const int TPZFGEOELEMENTRIANGLEID = 202;
/** @brief Id of geometric quadrilateral element */
const int TPZFGEOELEMENTQUADID = 203;
/** @brief Id of geometric cube element */
const int TPZFGEOELEMENTCUBEID = 204;
/** @brief Id of geometric prismal element */
const int TPZFGEOELEMENTPRISMID = 205;
/** @brief Id of geometric tetrahedral element */
const int TPZFGEOELEMENTTETRAID = 206;
/** @brief Id of geometric pyramidal element */
const int TPZFGEOELEMENTPYRAMID = 207;

/** @brief Id of geometric node */
const int TPZGEONODEID = 208;

/** @brief Id of geometric BC (boundary condition) element */
const int TPZGEOELBCID = 209;

/** @brief Id of geometric mesh */
const int TPZGEOMESHID = 210;

/** @brief Id of discontinuous computational element */
const int TPZCOMPELDISCID = 211;

/** @brief Id of agglomerated element (clusterred) */
const int TPZAGGLOMERATEELID = 212;

/** @brief Id of interface element */
const int TPZINTERFACEELEMENTID = 213;

/** @brief Id of point interpolated element */
const int TPZINTELPOINTID = 214;
/** @brief Id of linear interpolated element */
const int TPZINTELLINEARID = 215;
/** @brief Id of triangular interpolated element */
const int TPZINTELTRIANGLEID = 216;
/** @brief Id of quadrilateral interpolated element */
const int TPZINTELQUADID = 217;
/** @brief Id of cube interpolated element */
const int TPZINTELCUBEID = 218;
/** @brief Id of prismal interpolated element */
const int TPZINTELPRISMID = 219;
/** @brief Id of tetrahedral interpolated element */
const int TPZINTELTETRAID = 220;
/** @brief Id of pyramidal interpolated element */
const int TPZINTELPYRAMID = 221;

/** @brief Id of computational sub mesh */
const int TPZSUBCOMPMESHID = 222;
/** @brief Id of computational mesh */
const int TPZCOMPMESHID = 223;
/** @brief Id of computational mesh for flow */
const int TPZFLOWCOMPMESHID = 224;

/** @brief Id of geometric clone mesh */
const int TPZGEOCLONEMESHID = 225;
/** @brief Id of computational clone mesh */
const int TPZCOMPCLONEMESHID = 226;
/** @brief Id of computational clone mesh for parallelism */
const int TPZPARCOMPCLONEMESHID = 227;

/** @brief Id of geometric ref pattern point element */
const int TPZGEOELREFPATPOINTID = 228;
/** @brief Id of geometric ref pattern linear element */
const int TPZGEOELREFPATLINEARID = 229;
/** @brief Id of geometric ref pattern triangular element */
const int TPZGEOELREFPATTRIANGLEID = 230;
/** @brief Id of geometric ref pattern quadrilateral element */
const int TPZGEOELREFPATQUADID = 231;
/** @brief Id of geometric ref pattern cube element */
const int TPZGEOELREFPATCUBEID = 232;
/** @brief Id of geometric ref pattern prismal element */
const int TPZGEOELREFPATPRISMID = 233;
/** @brief Id of geometric ref pattern tetrahedral element */
const int TPZGEOELREFPATTETRAID = 234;
/** @brief Id of geometric ref pattern pyramidal element */
const int TPZGEOELREFPATPYRAMID = 235;
/** @brief Id of geometric ref pattern point mapped element */
const int TPZGEOELREFPATMAPPEDPOINTID = 236;
/** @brief Id of geometric ref pattern linear mapped element */
const int TPZGEOELREFPATMAPPEDLINEID = 237;
/** @brief Id of geometric ref pattern triangular mapped element */
const int TPZGEOELREFPATMAPPEDTRIANGLEID = 238;
/** @brief Id of geometric ref pattern quadrilateral mapped element */
const int TPZGEOELREFPATMAPPEDQUADRILATERALID = 239;
/** @brief Id of geometric ref pattern cube mapped element */
const int TPZGEOELREFPATMAPPEDCUBEID = 240;
/** @brief Id of geometric ref pattern prismal mapped element */
const int TPZGEOELREFPATMAPPEDPRISMID = 241;
/** @brief Id of geometric ref pattern tetrahedral mapped element */
const int TPZGEOELREFPATMAPPEDTETRAHEDRAID = 242;
/** @brief Id of geometric ref pattern pyramidal mapped element */
const int TPZGEOELREFPATMAPPEDPYRAMIDID = 243;

/** @brief Id of HDiv point element */
const int TPZHDIVPOINTID = 244;
/** @brief Id of HDiv linear element */
const int TPZHDIVLINEARID = 245;
/** @brief Id of HDiv triangular element */
const int TPZHDIVTRIANGLEID = 246;
/** @brief Id of HDiv quadrilateral element */
const int TPZHDIVQUADID = 247;
/** @brief Id of HDiv cube element */
const int TPZHDIVCUBEID = 248;
/** @brief Id of HDiv prismal element */
const int TPZHDIVPRISMID = 249;
/** @brief Id of HDiv tetrahedral element */
const int TPZHDIVTETRAID = 250;
/** @brief Id of HDiv pyramidal element */
const int TPZHDIVPYRAMID = 251;

/** @brief Id of HDiv boundary point element */
const int TPZHDIVBOUNDPOINTID = 252;
/** @brief Id of HDiv boundary linear element */
const int TPZHDIVBOUNDLINEARID = 253;
/** @brief Id of HDiv boundary triangular element */
const int TPZHDIVBOUNDTRIANGLEID = 254;
/** @brief Id of HDiv boundary quadrilateral element */
const int TPZHDIVBOUNDQUADID = 255;

/** @brief Id of HDiv two-dimensional boundary point element (?) */
const int TPZHDIVBOUND2POINTID = 256;
/** @brief Id of HDiv two-dimensional boundary linear element (?) */
const int TPZHDIVBOUND2LINEARID = 257;
/** @brief Id of HDiv two-dimensional boundary triangular element (?) */
const int TPZHDIVBOUND2TRIANGLEID = 258;
/** @brief Id of HDiv two-dimensional boundary quadrilateral element (?) */
const int TPZHDIVBOUND2QUADID = 259;

//-----
// /** @brief Id of multiphysics point element */
//const int TPZMULTIPHYSICSPOINTID = 260;
// /** @brief Id of multiphysics linear element */
//const int TPZMULTIPHYSICSLINEARID = 261;
// /** @brief Id of multiphysics triangular element */
//const int TPZMULTIPHYSICSTRIANGLEID = 262;
// /** @brief Id of multiphysics quadrilateral element */
//const int TPZMULTIPHYSICSQUADID = 263;
// /** @brief Id of multiphysics cube element */
//const int TPZMULTIPHYSICSCUBEID = 264;
// /** @brief Id of multiphysics prismal element */
//const int TPZMULTIPHYSICSPRISMID = 265;
// /** @brief Id of multiphysics tetrahedral element */
//const int TPZMULTIPHYSICSTETRAID = 266;
// /** @brief Id of multiphysics pyramidal element */
//const int TPZMULTIPHYSICSPYRAMID = 267;
//-------

/** @brief Id of geometric blended point element (?) */
const int TPZGEOBLENDPOINTID = 403;
/** @brief Id of geometric blended linear element (?) */
const int TPZGEOBLENDLINEARID = 404;
/** @brief Id of geometric blended quadrilateral element (?) */
const int TPZGEOBLENDQUADID = 405;
/** @brief Id of geometric blended triangular element (?) */
const int TPZGEOBLENDTRIANGLEID = 406;
/** @brief Id of geometric blended cube element (?) */
const int TPZGEOBLENDCUBEID = 407;
/** @brief Id of geometric blended prismal element (?) */
const int TPZGEOBLENDPRISMID = 408;
/** @brief Id of geometric blended pyramidal element (?) */
const int TPZGEOBLENDPYRAMIDID = 409;
/** @brief Id of geometric blended tetrahedral element (?) */
const int TPZGEOBLENDTETRAHEDRAID = 410;
/** @} */

/** @brief Function to register mesh classes */
void RegisterMeshClasses();
//template class TPZGeoElement<TPZShapePoint,TPZGeoPoint,TPZRefPoint>;
//template class TPZGeoElement<TPZShapeLinear,TPZGeoLinear,TPZRefLinear>;
//template class TPZGeoElement<TPZShapeTriang,TPZGeoTriangle,TPZRefTriangle>;
//template class TPZGeoElement<TPZShapeQuad,TPZGeoQuad,TPZRefQuad>;
//template class TPZGeoElement<TPZShapeCube,TPZGeoCube,TPZRefCube>;
//template class TPZGeoElement<TPZShapePrism,TPZGeoPrism,TPZRefPrism>;
//template class TPZGeoElement<TPZShapeTetra,TPZGeoTetrahedra,TPZRefTetrahedra>;
//template class TPZGeoElement<TPZShapePiram,TPZGeoPyramid,TPZRefPyramid>;


#endif
