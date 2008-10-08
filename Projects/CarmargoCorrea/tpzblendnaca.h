#ifndef TPZBLENDNACA_H
#define TPZBLENDNACA_H

#include <math.h>
#include "pzreal.h"
#include "pzvec.h"

/**
	@author caju2008 <caju@skol>
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

private:

    REAL P();
    REAL M();
    REAL TT();

    // Mean line for the wing
    REAL yc(REAL x);
    REAL dyc(REAL x);

    // thickness
    REAL yt(REAL x);

    // superior profile
    REAL xu(REAL x);
    REAL yu(REAL x);

    // inferior profile
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
