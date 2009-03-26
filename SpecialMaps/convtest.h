#ifndef CONVTEST_H
#define CONVTEST_H

#include <iostream>
#include <cstdlib>

#include "pzgmesh.h"
#include "pzvec.h"
#include "pzgnode.h"
#include "pzgeoel.h"
#include "pzadmchunk.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzmaterial.h"
#include "pzerror.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "TPZGeoElement.h"
#include "TPZRefPattern.h"
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
#include "tpzgeoelmapped.h"
#include "tpzcurvedtriangle.h"
#include "tpzquadratictrig.h"
#include "tpzquadratictetra.h"
#include "tpzarc3d.h"
#include "pzshapetriang.h"
#include "pzgeotriangle.h"
#include "pzlog.h"
#include "convtest.h"
#include <sstream>
#include "pzgeoelside.h"

using namespace std;
using namespace pzgeom;

  /**
  / Class made by Paulo Cesar de Alvarenga Lucci (Caju)
  / LabMeC - FEC - UNICAMP
  / 2007
 */

class ConvTest
{

public:

    ConvTest();
   ~ConvTest();

    void JacobianConv(TPZGeoEl &Object, TPZVec< REAL > QsiEta);
    void JacobianConv(TPZGeoElSide &Object, TPZVec< REAL > QsiEta);

};

#endif
