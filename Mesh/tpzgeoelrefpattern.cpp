/***************************************************************************
                          tpzgeoelrefpattern.cc  -  description
                             -------------------
    begin                : Tue Dec 23 2003
    copyright            : (C) 2003 by LabMeC - DES - FEC - UNICAMP (Edimar Cesar Rylo) & EMBRAER
    email                : cesar@labmec.fec.unicamp.br
 ***************************************************************************/

#include "tpzgeoelrefpattern.h"
#include "TPZGeoCube.h"
#include "pzshapecube.h"
#include "TPZRefCube.h"
#include "pzshapelinear.h"
#include "TPZGeoLinear.h"
#include "TPZRefLinear.h"
#include "pzrefquad.h"
#include "pzshapequad.h"
#include "pzgeoquad.h"
#include "pzshapetriang.h"
#include "pzreftriangle.h"
#include "pzgeotriangle.h"
#include "pzshapeprism.h"
#include "pzrefprism.h"
#include "pzgeoprism.h"
#include "pzshapetetra.h"
#include "pzreftetrahedra.h"
#include "pzgeotetrahedra.h"
#include "pzshapepiram.h"
#include "pzrefpyram.h"
#include "pzgeopyramid.h"
#include "pzrefpoint.h"
#include "pzgeopoint.h"
#include "pzshapepoint.h"
#include "pzgmesh.h"
#include "pzgeoel.h"
#include "TPZGeoElement.h"
#include "TPZRefPattern.h"
#include "pzvec.h"
#include "pzmanvector.h"

class TPZGeoElRefPattern<TPZShapeCube,TPZGeoCube>;
class TPZGeoElRefPattern<TPZShapeLinear,TPZGeoLinear>;
class TPZGeoElRefPattern<TPZShapeQuad,TPZGeoQuad>;
class TPZGeoElRefPattern<TPZShapeTriang,TPZGeoTriangle>;
class TPZGeoElRefPattern<TPZShapePrism,TPZGeoPrism>;
class TPZGeoElRefPattern<TPZShapeTetra,TPZGeoTetrahedra>;
class TPZGeoElRefPattern<TPZShapePiram,TPZGeoPyramid>;
class TPZGeoElRefPattern<TPZShapePoint,TPZGeoPoint>;
