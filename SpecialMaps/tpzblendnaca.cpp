/**
 * @file
 * @brief Contains the implementation of the TPZBlendNACA methods. 
 */
#include "tpzblendnaca.h"

#include "pzreal.h"
#include "pzvec.h"
#include "TPZGeoElement.h"
#include "pzshapequad.h"
#include "pzrefquad.h"
#include "TPZGeoLinear.h"

#include <iostream>
#include <fstream>

using namespace pzgeom;
using namespace pzshape;
using namespace pzrefine;

TPZBlendNACA::TPZBlendNACA()
{
}

TPZBlendNACA::TPZBlendNACA(REAL cord, int FourDigits, REAL angle, TPZVec<REAL> &x0) :
fCord(cord), fFourDigits(FourDigits), fAngle(angle)
{
    fX0[0] = x0[0];
    fX0[1] = x0[1];
    fX0[2] = x0[2];
    fP = P();
    fM = M();
    fTT = TT();
}

TPZBlendNACA::~TPZBlendNACA()
{
}

REAL TPZBlendNACA::P()
{
    int aux = fFourDigits/100;
    aux -= ((int)(aux/10))*10;
    return (REAL)(aux/10.);
}

REAL TPZBlendNACA::M()
{
    int aux = fFourDigits/1000;
    return (REAL)(aux/100.)*fCord;
}

REAL TPZBlendNACA::TT()
{
    int aux = fFourDigits - ((int)(fFourDigits/100))*100;
    return (REAL)(aux/100.)*fCord;
}

REAL TPZBlendNACA::yc(REAL x)
{
    if(x/fCord<fP)
    {
        return fM/fP/fP*(2.*fP*x/fCord-x*x/fCord/fCord);
    }
    else
    {
        return fM/(1.-fP)/(1.-fP)*(1.-2.*fP+2.*fP*x/fCord-x*x/fCord/fCord);
    }
}

REAL TPZBlendNACA::dyc(REAL x)
{
    if(x/fCord<fP)
    {
        return 2.*fM/fP/fP*(fP-x/fCord)/fCord;
    }
    else
    {
        return 2.*fM/(1.-fP)/(1.-fP)*(fP-x/fCord)/fCord;
    }
}

REAL TPZBlendNACA::yt(REAL x)
{
    REAL aux = x/fCord;
    const REAL a0 =  1.4845,
	a1 = -0.6300,
	a2 = -1.7580,
	a3 =  1.4215,
	a4 = -0.5075;
	
    return fTT * (a0*sqrt(aux) + a1*aux + a2*aux*aux + a3*aux*aux*aux + a4*aux*aux*aux*aux);
}

REAL TPZBlendNACA::xu(REAL x)
{
    return x-yt(x)*sin(atan(dyc(x)));
}

REAL TPZBlendNACA::yu(REAL x)
{
    return yc(x) + yt(x)*cos(atan(dyc(x)));
}

REAL TPZBlendNACA::xl(REAL x)
{
    return x+yt(x)*sin(atan(dyc(x)));
}

REAL TPZBlendNACA::yl(REAL x)
{
    return yc(x) - yt(x)*cos(atan(dyc(x)));
}

REAL TPZBlendNACA::NearestParameter(TPZVec<REAL> &pt, int &uplow, int maxPt)
{
    REAL distminlow = 20.*fCord;
    REAL distminup  = 20.*fCord;
    REAL distlow,distup;
    REAL ptl[2],ptu[2];
    int ip,maxp=maxPt;
    REAL par,parlow;
	REAL parup = 0.0L;
    for(ip=0; ip<=maxp; ip++)
    {
        par = ip*fCord/maxp;
        ptu[0] = xua(par);
        ptu[1] = yua(par);
        ptl[0] = xla(par);
        ptl[1] = yla(par);
        distlow = (ptl[0]-pt[0])*(ptl[0]-pt[0])+(ptl[1]-pt[1])*(ptl[1]-pt[1]);
        distup  = (ptu[0]-pt[0])*(ptu[0]-pt[0])+(ptu[1]-pt[1])*(ptu[1]-pt[1]);
        if(distlow < distminlow) parlow = par;
        if(distup < distminup) parup = par;
        distminlow = distminlow < distlow ? distminlow : distlow;
        distminup = distminup < distup ? distminup : distup;
    }
    REAL delpar = (0.1L)/maxp;
    if(distminlow < distminup)
    {
        uplow = 0;
        REAL distprev = distminlow;
        par = parlow;
        while(fabs(delpar) > 0.00001/maxp)
        {
            ptl[0] = xla(par+delpar);
            ptl[1] = yla(par+delpar);
            distlow = (ptl[0]-pt[0])*(ptl[0]-pt[0])+(ptl[1]-pt[1])*(ptl[1]-pt[1]);
            if(distlow < distprev)
            {
                par += delpar;
                distprev = distlow;
            }
            else if (delpar > 0.)
            {
                delpar *= -1.;
            }
            else
            {
                delpar *= -0.1;
            }
        }
    } 
    else 
    {
        uplow = 1;
        REAL distprev = distminup;
        par = parup;
        while(fabs(delpar) > 0.001/maxp)
        {
            ptu[0] = xua(par+delpar);
            ptu[1] = yua(par+delpar);
            distup = (ptu[0]-pt[0])*(ptu[0]-pt[0])+(ptu[1]-pt[1])*(ptu[1]-pt[1]);
            if(distup < distprev)
            {
                par += delpar;
                distprev = distup;
            }
            else if (delpar > 0.)
            {
                delpar *= -1.;
            }
            else
            {
                delpar *= -0.1;
            }
        }
    }
    return par;
}

///************************************
REAL TPZBlendNACA::xua(REAL x)
{
    return fX0[0]+(xu(x)-fCord/2.)*cos(fAngle) + yu(x) * sin(fAngle) + fCord/2.;
}

REAL TPZBlendNACA::yua(REAL x)
{
    return fX0[1]+yu(x)*cos(fAngle) - (xu(x)-fCord/2.) * sin(fAngle);
}

REAL TPZBlendNACA::xla(REAL x)
{
    return fX0[0]+(xl(x)-fCord/2.)*cos(fAngle) + yl(x) * sin(fAngle) + fCord/2.;
}

REAL TPZBlendNACA::yla(REAL x)
{
    return fX0[1]+yl(x)*cos(fAngle) - (xl(x)-fCord/2.) * sin(fAngle);
}
///************************************
void TPZBlendNACA::ProjectPoint(TPZVec<REAL> &pt, int maxPt)
{
    int uplow;
    REAL par = NearestParameter(pt,uplow, maxPt);
    if(uplow == 0)
    {
        pt[0] = xla(par);
        pt[1] = yla(par);
    }
    else
    {
        pt[0] = xua(par);
        pt[1] = yua(par);
    }
}

#include "tpzgeoelmapped.h"

/**
 * Creates a geometric element according to the type of the father element
 */

TPZGeoEl *TPZBlendNACA::CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
										 TPZVec<int>& nodeindexes,
										 int matid,
										 int& index)
{
	return CreateGeoElementMapped(mesh,type,nodeindexes,matid,index);
}


