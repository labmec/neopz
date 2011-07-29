#ifndef TPZCHANGEEL_H
#define TPZCHANGEEL_H

#include "pzgeotriangle.h"
#include "tpzgeoblend.h"
#include "pzgeoel.h"
#include "pznoderep.h"

#include <iostream>

  /**
  / Class made by Paulo Cesar de Alvarenga Lucci (Caju)
  / LabMeC - FEC - UNICAMP
  / 2007
   */

/**
 * @author Paulo Cesar de Alvarenga Lucci (Caju)
 * @ingroup geometry
 */
class TPZChangeEl {
public:

    TPZChangeEl();
   ~TPZChangeEl();

    /** @brief Turns an linear triangle or quadrilateral to quadratic */
    static TPZGeoEl * ChangeToQuadratic(TPZGeoMesh *Mesh, int ElemIndex);

    /** @brief Turns a regular element into a geoblend */
    static TPZGeoEl * ChangeToGeoBlend(TPZGeoMesh *Mesh, int ElemIndex);

    /** @brief Turns an quadratic triangle or quadrilateral to linear */
    static TPZGeoEl * ChangeToLinear(TPZGeoMesh *Mesh, int ElemIndex);

    /** @brief Slide correct nodes of an triangle or quadrilateral to the quarterpoint with respect to a given side */
    static TPZGeoEl * QuarterPoints(TPZGeoMesh *Mesh, int ElemIndex, int side);

private:
    static void AdjustNeighbourhood(TPZGeoEl* OldElem, TPZGeoEl*NewElem);

};

#endif
