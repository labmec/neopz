//
// C++ Implementation: tpzgeoelmapped
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "tpzgeoelmapped.h"

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
#include "tpzgeoelrefpattern.h"

using namespace pzgeom;
using namespace pzshape;

template class TPZGeoElMapped< TPZGeoElRefPattern<TPZGeoCube> >;
template class TPZGeoElMapped< TPZGeoElement<TPZGeoCube,pzrefine::TPZRefCube> >;

template<>
int TPZGeoElMapped< TPZGeoElRefPattern<TPZGeoCube> >::ClassId() const {
  return 200;
    }

    template<>
        int TPZGeoElMapped< TPZGeoElement<TPZGeoCube,pzrefine::TPZRefCube> >::ClassId() const {
      return 200;
        }
