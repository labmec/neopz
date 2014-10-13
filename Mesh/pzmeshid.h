/**
 * @file
 * @brief Contains the Ids for several geometric structures.
 */

#ifndef PZMESHIDH
#define PZMESHIDH

/** \addtogroup geometry
 * @{
 */


/** @brief Identifier indicating the no material is associated */
const int GMESHNOMATERIAL = -9999;

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
const int TPZHDIVPRESSUREBOUNDPOINTID = 252;
/** @brief Id of HDiv boundary linear element */
const int TPZHDIVPRESSUREBOUNDLINEARID = 253;
/** @brief Id of HDiv boundary triangular element */
const int TPZHDIVPRESSUREBOUNDTRIANGLEID = 254;
/** @brief Id of HDiv boundary quadrilateral element */
const int TPZHDIVPRESSUREBOUNDQUADID = 255;

/** @brief Id of HDiv point element */
const int TPZHDIVPOINT_PRESSURE_ID = 260;
/** @brief Id of HDiv linear element */
const int TPZHDIVLINEAR_PRESSURE_ID = 261;
/** @brief Id of HDiv triangular element */
const int TPZHDIVTRIANGLE_PRESSURE_ID = 262;
/** @brief Id of HDiv quadrilateral element */
const int TPZHDIVQUAD_PRESSURE_ID = 263;
/** @brief Id of HDiv cube element */
const int TPZHDIVCUBE_PRESSURE_ID = 264;
/** @brief Id of HDiv prismal element */
const int TPZHDIVPRISM_PRESSURE_ID = 265;
/** @brief Id of HDiv tetrahedral element */
const int TPZHDIVTETRA_PRESSURE_ID = 266;
/** @brief Id of HDiv pyramidal element */
const int TPZHDIVPYRAM_PRESSURE_ID = 267;

/** @brief Id of HDivFull point element */
const int TPZHDIVPOINT_FULL_ID = 268;
/** @brief Id of HDivFull linear element */
const int TPZHDIVLINEAR_FULL_ID = 269;
/** @brief Id of HDivFull triangular element */
const int TPZHDIVTRIANGLE_FULL_ID = 270;
/** @brief Id of HDivFull quadrilateral element */
const int TPZHDIVQUAD_FULL_ID = 271;
/** @brief Id of HDivFull cube element */
const int TPZHDIVCUBE_FULL_ID = 272;
/** @brief Id of HDivFull prismal element */
const int TPZHDIVPRISM_FULL_ID = 273;
/** @brief Id of HDivFull tetrahedral element */
const int TPZHDIVTETRA_FULL_ID = 274;
/** @brief Id of HDivFull pyramidal element */
const int TPZHDIVPYRAM_FULL_ID = 275;



/** @brief Id of HDiv one-dimensional boundary point element */
const int TPZHDIVBOUND2POINTID = 256;
/** @brief Id of HDiv two-dimensional boundary linear element  */
const int TPZHDIVBOUND2LINEARID = 257;
/** @brief Id of HDiv three-dimensional boundary triangular element */
const int TPZHDIVBOUND2TRIANGLEID = 258;
/** @brief Id of HDiv three-dimensional boundary quadrilateral element */
const int TPZHDIVBOUND2QUADID = 259;

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

/** @brief Id of computational memory point element (?) */
const int TPZCOMPELWITHMEMPOINTID = 411;
/** @brief Id of computational memory linear element (?) */
const int TPZCOMPELWITHMEMLINEARID = 412;
/** @brief Id of computational memory triangle element (?) */
const int TPZCOMPELWITHMEMTRIANGID = 413;
/** @brief Id of computational memory quadrilateral element (?) */
const int TPZCOMPELWITHMEMQUADID = 414;
/** @brief Id of computational memory cube element (?) */
const int TPZCOMPELWITHMEMCUBEID = 415;
/** @brief Id of computational memory tetrahedral element (?) */
const int TPZCOMPELWITHMEMTETRAID = 416;
/** @brief Id of computational memory prismal element (?) */
const int TPZCOMPELWITHMEMPRISMID = 417;
/** @brief Id of computational memory pyramidal element (?) */
const int TPZCOMPELWITHMEMPIRAMID = 418;

const int TPZCOMPMESHREFERREDID = 419;

//For Multiphysics
/** @brief Id of computational memory point element (?) */
const int TPZMPCOMPELWITHMEMPOINTID = 420;
/** @brief Id of computational memory linear element (?) */
const int TPZMPCOMPELWITHMEMLINEARID = 421;
/** @brief Id of computational memory triangle element (?) */
const int TPZMPCOMPELWITHMEMTRIANGID = 422;
/** @brief Id of computational memory quadrilateral element (?) */
const int TPZMPCOMPELWITHMEMQUADID = 423;
/** @brief Id of computational memory cube element (?) */
const int TPZMPCOMPELWITHMEMCUBEID = 424;
/** @brief Id of computational memory tetrahedral element (?) */
const int TPZMPCOMPELWITHMEMTETRAID = 425;
/** @brief Id of computational memory prismal element (?) */
const int TPZMPCOMPELWITHMEMPRISMID = 426;
/** @brief Id of computational memory pyramidal element (?) */
const int TPZMPCOMPELWITHMEMPIRAMID = 427;

//For MultiPhysicsInterface
/** @brief Id of computational memory interface multiphysics element (?) */
const int TPZMPCOMPELWITHMEMINTERFACE = 428;

/** @} */

#endif
