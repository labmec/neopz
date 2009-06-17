// -*- c++ -*-
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

const int TPZFGEOELEMENTPOINTID = 200;
const int TPZFGEOELEMENTLINEARID = 201;
const int TPZFGEOELEMENTRIANGLEID = 202;
const int TPZFGEOELEMENTQUADID = 203;
const int TPZFGEOELEMENTCUBEID = 204;
const int TPZFGEOELEMENTPRISMID = 205;
const int TPZFGEOELEMENTTETRAID = 206;
const int TPZFGEOELEMENTPYRAMID = 207;

const int TPZGEONODEID = 208;

const int TPZGEOELBCID = 209;

const int TPZGEOMESHID = 210;

const int TPZCOMPELDISCID = 211;

const int TPZAGGLOMERATEELID = 212;

const int TPZINTERFACEELEMENTID = 213;

const int TPZINTELPOINTID = 214;
const int TPZINTELLINEARID = 215;
const int TPZINTELTRIANGLEID = 216;
const int TPZINTELQUADID = 217;
const int TPZINTELCUBEID = 218;
const int TPZINTELPRISMID = 219;
const int TPZINTELTETRAID = 220;
const int TPZINTELPYRAMID = 221;

const int TPZSUBCOMPMESHID = 222;
const int TPZCOMPMESHID = 223;
const int TPZFLOWCOMPMESHID = 224;

const int TPZGEOCLONEMESHID = 225;
const int TPZCOMPCLONEMESHID = 226;
const int TPZPARCOMPCLONEMESHID = 227;

const int TPZGEOELREFPATPOINTID = 228;
const int TPZGEOELREFPATLINEARID = 229;
const int TPZGEOELREFPATTRIANGLEID = 230;
const int TPZGEOELREFPATQUADID = 231;
const int TPZGEOELREFPATCUBEID = 232;
const int TPZGEOELREFPATPRISMID = 233;
const int TPZGEOELREFPATTETRAID = 234;
const int TPZGEOELREFPATPYRAMID = 235;
const int TPZGEOELREFPATMAPPEDPOINTID = 236;
const int TPZGEOELREFPATMAPPEDLINEID = 237;
const int TPZGEOELREFPATMAPPEDTRIANGLEID = 238;
const int TPZGEOELREFPATMAPPEDQUADRILATERALID = 239;
const int TPZGEOELREFPATMAPPEDCUBEID = 240;
const int TPZGEOELREFPATMAPPEDPRISMID = 241;
const int TPZGEOELREFPATMAPPEDTETRAHEDRAID = 242;
const int TPZGEOELREFPATMAPPEDPYRAMIDID = 243;

const int TPZHDIVPOINTID = 244;
const int TPZHDIVLINEARID = 245;
const int TPZHDIVTRIANGLEID = 246;
const int TPZHDIVQUADID = 247;
const int TPZHDIVCUBEID = 248;
const int TPZHDIVPRISMID = 249;
const int TPZHDIVTETRAID = 250;
const int TPZHDIVPYRAMID = 251;

const int TPZGEOBLENDPOINTID = 403;
const int TPZGEOBLENDLINEARID = 404;
const int TPZGEOBLENDQUADID = 405;
const int TPZGEOBLENDTRIANGLEID = 406;
const int TPZGEOBLENDCUBEID = 407;
const int TPZGEOBLENDPRISMID = 408;
const int TPZGEOBLENDPYRAMIDID = 409;
const int TPZGEOBLENDTETRAHEDRAID = 410;

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
