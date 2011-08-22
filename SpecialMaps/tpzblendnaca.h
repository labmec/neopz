/**
 * @file
 * @brief Contains the TPZBlendNACA class. It is a special map.
 */
#ifndef TPZBLENDNACA_H
#define TPZBLENDNACA_H

#include <math.h>
#include "pzreal.h"
#include "pzvec.h"
#include "pzeltype.h"

class TPZGeoEl;
class TPZGeoMesh;
/**
 * @author caju2008 < caju\@skol >
 * @ingroup geometry
 * @brief Special map to NACA. \ref geometry "Geometry"
 */
class TPZBlendNACA
{
public:
	
    TPZBlendNACA();
    TPZBlendNACA(REAL cord, int FourDigits, REAL angle, TPZVec<REAL> &x0);
	~TPZBlendNACA();
	
    /// with attack angle
    /// superior profile
    REAL xua(REAL x);
    REAL yua(REAL x);
	
    /// inferior profile
    REAL xla(REAL x);
    REAL yla(REAL x);
	
    void ProjectPoint(TPZVec<REAL> &pt, int maxPt = 1000);
	
public:
	/* * @brief Creates a geometric element according to the type of the father element */
	static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
									  TPZVec<int>& nodeindexes,
									  int matid,
									  int& index);

private:
	
    REAL P();
    REAL M();
    REAL TT();
	
    /** @brief Mean line for the wing */
    REAL yc(REAL x);
    REAL dyc(REAL x);
	
    /** @brief Thickness */
    REAL yt(REAL x);
	
    /** @brief Superior profile */
    REAL xu(REAL x);
    REAL yu(REAL x);
	
    /** @brief Inferior profile */
    REAL xl(REAL x);
    REAL yl(REAL x);
	
    REAL NearestParameter(TPZVec<REAL> &pt, int &uplow, int maxPt);
	
public:
    REAL fCord;
    int  fFourDigits;
    REAL fAngle;
    REAL fX0[3];
    REAL fP;
    REAL fM;
    REAL fTT;
};

#endif
