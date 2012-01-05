/**
 * @file
 * @brief Contains the TPZChangeEl class. It is a special map.
 */
#ifndef TPZCHANGEEL_H
#define TPZCHANGEEL_H

#include "pzgeotriangle.h"
#include "tpzgeoblend.h"
#include "pzgeoel.h"
#include "pznoderep.h"

#include <iostream>

/**
 * @author Paulo Cesar de Alvarenga Lucci (Caju)
 * @ingroup geometry
 * @since 2007
 * @brief Special map .... \ref geometry "Geometry"
 */
class TPZChangeEl {
public:

    TPZChangeEl();
   ~TPZChangeEl();

    /** @brief Turns an linear geoelement to quadratic */
    static TPZGeoEl * ChangeToQuadratic(TPZGeoMesh *Mesh, int ElemIndex);

    /** @brief Turns a regular element into a geoblend */
    static TPZGeoEl * ChangeToGeoBlend(TPZGeoMesh *Mesh, int ElemIndex);

    /** @brief Slide middle nodes of an quadratic geoelement to the quarterpoint with respect to a given side */
    static TPZGeoEl * ChangeToQuarterPoint(TPZGeoMesh *Mesh, int ElemIndex, int targetSide);
};

#endif
