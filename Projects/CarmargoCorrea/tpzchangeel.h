#ifndef TPZCHANGEEL_H
#define TPZCHANGEEL_H

#include "pzgeotriangle.h"
#include "tpzgeoblend.h"
#include "pzgeoel.h"
#include "pznoderep.h"

#include <iostream>

using namespace std;
using namespace pzgeom;
using namespace pztopology;

  /**
  / Class made by Paulo Cesar de Alvarenga Lucci (Caju)
  / LabMeC - FEC - UNICAMP
  / 2007
 */

class TPZChangeEl {
public:

    TPZChangeEl();
   ~TPZChangeEl();

    /** Turns an linear triangle or quadrilateral to quadratic */
    void ChangeToQuadratic(TPZGeoMesh *Mesh, int ElemIndex);

    /** Turns an quadratic triangle or quadrilateral to linear */
    void ChangeToLinear(TPZGeoMesh *Mesh, int ElemIndex);

    /** Slide correct nodes of an triangle or quadrilateral to the quarterpoint with respect to a given side */
    void QuarterPoints(TPZGeoMesh *Mesh, int ElemIndex, int side);

};

#endif
