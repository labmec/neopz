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
 * @brief Special map. It is util to convert a linear element for quadratic element, but the same topology. \ref geometry "Geometry"
 */
/** Also turns slide middle nodes of an quadratic element to the quarterpoint with respect to a indicated side. */
class TPZChangeEl {
public:

    TPZChangeEl();
   ~TPZChangeEl();

    /** @brief Turns an linear geoelement to quadratic */
    static TPZGeoEl * ChangeToQuadratic(TPZGeoMesh *Mesh, int64_t ElemIndex);

    /** @brief Turns a regular element into a geoblend */
    static TPZGeoEl * ChangeToGeoBlend(TPZGeoMesh *Mesh, int64_t ElemIndex);

    /** @brief Slide middle nodes of an quadratic geoelement to the quarterpoint with respect to a given side */
    static TPZGeoEl * ChangeToQuarterPoint(TPZGeoMesh *Mesh, int64_t ElemIndex, int targetSide);

    /**
	 * @brief Return if a given point x is near to some node of a given geo element
	 * @param gel [in] given geo element
	 * @param x [in] given point
	 * @param meshNode [out] id of node that is in the x range
	 * @param tol [in] x range radius
	 */
	static bool NearestNode(TPZGeoEl * gel, TPZVec<REAL> &x, int64_t &meshNode, double tol);
    /**
	 * @brief Return the id of the node into the geometric mesh nearest a given point x
	 * @param gmesh [in] given geometric mesh
	 * @param x [in] given point
	 * @param tol [in] x range radius
	 */
    static int64_t NearestNode(TPZGeoMesh * gmesh, TPZVec<REAL> &x, double tol);
    static bool CreateMiddleNodeAtEdge(TPZGeoMesh *Mesh, int64_t ElemIndex, int edge, int64_t &middleNodeId);
};

#endif
