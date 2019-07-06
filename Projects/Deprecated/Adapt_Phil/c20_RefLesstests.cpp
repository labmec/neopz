/***************************************************************************
                          c20_RefLesstests.cpp  -  description
                             -------------------
    begin                : Fri Dec 5 2003
    copyright            : (C) 2003 by cesar
    email                : cesar@becks.fec.unicamp.br
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

REAL Coord[24][3] = {
  {0.,0.,0.},{1.,0.,0.},{2.,0.,0.},{3.,0.,0.}, //y=0,z=0
  {0.,1.,0.},{1.,1.,0.},{2.,1.,0.},{3.,1.,0.}, //y=1,z=0

  {0.,0.,1.},{1.,0.,1.},{2.,0.,1.},{3.,0.,1.}, //y=0,z=1
  {0.,1.,1.},{1.,1.,1.},{2.,1.,1.},{3.,1.,1.}, //y=1,z=1

  {0.,0.,2.},{1.,0.,2.},{2.,0.,2.},{3.,0.,2.}, //y=0,z=2
  {0.,1.,2.},{1.,1.,2.},{2.,1.,2.},{3.,1.,2.} //y=1,z=2
};

int NodePerEl [14] = {8,8,5,5,5,6,6,4,4,4,4,4,6,6};
int Connects [] = {
  0,1,5,4,
  8,9,13,12,
  1,2,6,5,14,
  1,2,10,9,14,
  1,5,13,9,14,
  

  
int newteste(){
  TPZGeoElRefLess <TPZShapeCube,TPZGeoCube> el1;
  TPZGeoElRefLess <TPZShapeLinear,TPZGeoLinear> el2;
  TPZGeoElRefLess <TPZShapeQuad,TPZGeoQuad> el3;
  TPZGeoElRefLess <TPZShapeTriang,TPZGeoTriangle> el4;
  TPZGeoElRefLess <TPZShapePrism,TPZGeoPrism> el5;
  TPZGeoElRefLess <TPZShapeTetra,TPZGeoTetrahedra> el6;
  TPZGeoElRefLess <TPZShapePiram,TPZGeoPyramid> el7;
  TPZGeoElRefLess <TPZShapePoint,TPZGeoPoint> el8;
return 0;
}
