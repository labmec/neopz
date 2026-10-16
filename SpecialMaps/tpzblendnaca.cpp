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
#include "tpzgeoelrefpattern.h"

#include <iostream>
#include <fstream>

using namespace pzgeom;
using namespace pzshape;
using namespace pzrefine;

TPZBlendNACA::TPZBlendNACA() : fCord(1.), fFourDigits(0), fAngle(0.)
{
    fX0[0] = 0.;
    fX0[1] = 0.;
    fX0[2] = 0.;
}

TPZBlendNACA::TPZBlendNACA(REAL cord, int FourDigits, REAL angle, TPZVec<REAL> &x0) :
fCord(cord), fFourDigits(FourDigits), fAngle(angle)
{
    fX0[0] = x0[0];
    fX0[1] = x0[1];
    fX0[2] = x0[2];
    fP = P<REAL>();
    fM = M<REAL>();
    fTT = TT<REAL>();
}

TPZBlendNACA::~TPZBlendNACA()
{
}

template<class Type>
Type TPZBlendNACA::P()
{
    int aux = fFourDigits/100;
    aux -= ((int)(aux/10))*10;
    return (Type)(aux/10.);
}
template REAL TPZBlendNACA::P();

template<class Type>
Type TPZBlendNACA::M()
{
    int aux = fFourDigits/1000;
    return (Type)(aux/100.)*fCord;
}
template REAL TPZBlendNACA::M();


template<class toto>
toto TPZBlendNACA::TT()
{
    int aux = fFourDigits - ((int)(fFourDigits/100))*100;
    toto result;
    result = (toto)(aux/100.)*fCord;
    // std::cout << "TT = " << result << std::endl;
    return result;
}
template REAL TPZBlendNACA::TT();


template <class toto>
toto TPZBlendNACA::yc(toto x) const
{
    if(shapeFAD::val(x)/fCord<fP)
    {
        return fM/fP/fP*(2.*fP*x/fCord-x*x/fCord/fCord);
    }
    else
    {
        return fM/(1.-fP)/(1.-fP)*(1.-2.*fP+2.*fP*x/fCord-x*x/fCord/fCord);
    }
}

template REAL TPZBlendNACA::yc(REAL x) const;
template Fad<REAL> TPZBlendNACA::yc(Fad<REAL> x) const;

//*******
template <class toto>
toto TPZBlendNACA::dyc(toto x) const
{
    if(shapeFAD::val(x)/fCord<fP)
    {
        return 2.*fM/fP/fP*(fP-x/fCord)/fCord;
    }
    else
    {
        return 2.*fM/(1.-fP)/(1.-fP)*(fP-x/fCord)/fCord;
    }
}
template REAL TPZBlendNACA::dyc(REAL x) const;
template Fad<REAL> TPZBlendNACA::dyc(Fad<REAL> x) const;

//*******
template <class toto>
toto TPZBlendNACA::dtgphi(toto x) const
{
    if(shapeFAD::val(x)/fCord<fP)
    {
        return 2.*fM/fP/fP*(-1./fCord)/fCord;
    }
    else
    {
        return 2.*fM/(1.-fP)/(1.-fP)*(-1./fCord)/fCord;
    }
}
template REAL TPZBlendNACA::dtgphi(REAL x) const;
template Fad<REAL> TPZBlendNACA::dtgphi(Fad<REAL> x) const;

//********
template <class toto>
toto TPZBlendNACA::yt(toto x) const
{
    toto aux = x/fCord;
    const REAL a0 =  1.4845,
	a1 = -0.6300,
	a2 = -1.7580,
	a3 =  1.4215,
//	a4 = -0.5075;
	a4 = -0.518;
	
    toto val = fTT * (a0*sqrt(aux) + a1*aux + a2*aux*aux + a3*aux*aux*aux + a4*aux*aux*aux*aux);
    //std::cout << " x " << x << " yt = " << val << std::endl;
    return val;
}

//********
template <class toto>
toto TPZBlendNACA::dyt(toto x) const
{
    toto aux = x/fCord;
    const REAL a0 =  1.4845,
	a1 = -0.6300,
	a2 = -1.7580,
	a3 =  1.4215,
//	a4 = -0.5075;
	a4 = -0.518;
	
    return fTT/fCord * (a0/2./sqrt(aux) + a1 + 2.*a2*aux + 3.*a3*aux*aux + 4.*a4*aux*aux*aux);
}

template REAL TPZBlendNACA::dyt(REAL x) const;

template <class toto>
toto TPZBlendNACA::xl(toto x) const
{
    return x+yt(x)*sin(atan(tgphi(x)));
}
template REAL TPZBlendNACA::xl(REAL x) const;


//********
template <class toto>
toto TPZBlendNACA::xu(toto x) const
{
    toto ytval = yt(x);
    toto tgphival = tgphi(x);
    toto val = x-ytval*sin(atan(tgphival));
    // std::cout << "xu x " << x << " yt " << ytval << " tgphi " << tgphival << " val = " << val << std::endl;
    return val;
}

template REAL TPZBlendNACA::xu(REAL x) const;
template Fad<REAL> TPZBlendNACA::xu(Fad<REAL> x) const;
//*********

//********
template <class toto>
toto TPZBlendNACA::dxu(toto x) const
{
    auto loctgphi = tgphi(x);
    auto ytval = yt(x);
    auto dtgphival = dtgphi(x);
    auto dytval = dyt(x);
    auto resp =  1.-dytval*sin(atan(loctgphi))
        -ytval/(1+loctgphi*loctgphi)/sqrt(1.+loctgphi*loctgphi)*dtgphival;
    // std::cout << "x " << x << " tgphi " << loctgphi << " yt " << ytval << " dtgphi " << dtgphival << " dyt " << dytval << " resp = " << resp << std::endl;
    return resp;
}

template REAL TPZBlendNACA::dxu(REAL x) const;
template Fad<REAL> TPZBlendNACA::dxu(Fad<REAL> x) const;
//*********


template <class toto>
toto TPZBlendNACA::yu(toto x) const
{
    return yc(x) + yt(x)*cos(atan(tgphi(x)));
}
template REAL TPZBlendNACA::yu(REAL x) const;
//********

template <class toto>
toto TPZBlendNACA::dyu(toto x) const
{
    auto loctgphi = tgphi(x);
    auto resp =  loctgphi+dyt(x)*cos(atan(dyc(x)))
        -yt(x)*loctgphi/(1+loctgphi*loctgphi)/sqrt(1.+loctgphi*loctgphi)*dtgphi(x);
    return resp;
}
template REAL TPZBlendNACA::dyu(REAL x) const;
//********

template <class toto>
toto TPZBlendNACA::dxl(toto x) const
{
    toto loctgphi = tgphi(x);
    return 1.+dyt(x)*sin(atan(loctgphi))+
        yt(x)/(1+loctgphi*loctgphi)/sqrt(1.+loctgphi*loctgphi)*dtgphi(x);
}
template REAL TPZBlendNACA::dxl(REAL x) const;

//********
template <class toto>
toto TPZBlendNACA::yl(toto x) const
{
    return yc(x) - yt(x)*cos(atan(tgphi(x)));
}
template REAL TPZBlendNACA::yl(REAL x) const;

//********
template <class toto>
toto TPZBlendNACA::dyl(toto x) const
{
    toto loctgphi = tgphi(x);
    return loctgphi - dyt(x)*cos(atan(loctgphi))+
        yt(x)*loctgphi/(1+loctgphi*loctgphi)/sqrt(1.+loctgphi*loctgphi)*dtgphi(x);
}
template REAL TPZBlendNACA::dyl(REAL x) const;


//********

template<class toto>
void TPZBlendNACA::NearestParameter(TPZVec<REAL> &pt, int &uplow, int maxPt,toto& Result)
{
    REAL distminlow = 20.*fCord;
    REAL distminup  = 20.*fCord;
    REAL distlow,distup;
    REAL ptl[2],ptu[2];
    int ip,maxp=maxPt;
    REAL par,parlow = (REAL)0.0;
	REAL parup = (REAL)0.0;
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
    Result = par;
}
template<> void TPZBlendNACA::NearestParameter(TPZVec<REAL> &pt, int &uplow, int maxPt,REAL& Result);

///************************************ public methods
template <class toto>
toto TPZBlendNACA::xua(toto x) const
{
    return fX0[0]+(xu(x)-fCord/2.)*cos(fAngle) + yua(x) * sin(fAngle) + fCord/2.;
}

template REAL TPZBlendNACA::xua(REAL x) const;

//***************

template <class toto>
toto TPZBlendNACA::yua(toto x) const
{
    return fX0[1]+yu(x)*cos(fAngle) - (xu(x)-fCord/2.) * sin(fAngle);
}
template REAL TPZBlendNACA::yua(REAL x) const;

//************
template <class toto>
toto TPZBlendNACA::xla(toto x) const
{
    return fX0[0]+(xl(x)-fCord/2.)*cos(fAngle) + yl(x) * sin(fAngle) + fCord/2.;
}
template REAL TPZBlendNACA::xla(REAL x) const;

//************
template <class toto>
toto TPZBlendNACA::yla(toto x) const
{
    return fX0[1]+yl(x)*cos(fAngle) - (xl(x)-fCord/2.) * sin(fAngle);
}
template REAL TPZBlendNACA::yla(REAL) const;

template <class toto>
toto TPZBlendNACA::dxua(toto x) const
{
    return fX0[0]+(dxu(x))*cos(fAngle) + dyu(x) * sin(fAngle);
}
template REAL TPZBlendNACA::dxua(REAL x) const;

//***************

template <class toto>
toto TPZBlendNACA::dyua(toto x) const
{
    return dyu(x)*cos(fAngle) - (dxu(x)) * sin(fAngle);
}
template REAL TPZBlendNACA::dyua(REAL x) const;

//************
template <class toto>
toto TPZBlendNACA::dxla(toto x) const
{
    return (dxl(x))*cos(fAngle) + dyl(x) * sin(fAngle);
}
template REAL TPZBlendNACA::dxla(REAL x) const;

//************
template <class toto>
toto TPZBlendNACA::dyla(toto x) const
{
    return dyl(x)*cos(fAngle) - (xl(x)) * sin(fAngle);
}
template REAL TPZBlendNACA::dyla(REAL) const;
//template REAL TPZBlendNACA::yla(REAL x);

///************************************
template <class toto>
void TPZBlendNACA::ProjectPoint(TPZVec<toto> &pt, int maxPt)
{
    int uplow;
    
    toto par;
    NearestParameter(pt,uplow, maxPt,par);
    
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
template<> void TPZBlendNACA::ProjectPoint(TPZVec<REAL> &pt, int maxPt);

int TPZBlendNACA::ClassId() const{
    return Hash("TPZBlendNACA") ^ pzgeom::TPZNodeRep<2,pztopology::TPZLine>::ClassId() << 1;
}

#include "tpzgeoelmapped.h"

/**
 * Creates a geometric element according to the toto of the father element
 */

// TPZGeoEl *TPZBlendNACA::CreateGeoElement(TPZGeoMesh &mesh, MElementtoto toto,
// 										 TPZVec<int64_t>& nodeindexes,
// 										 int matid,
// 										 int64_t& index)
// {
// 	return CreateGeoElementMapped(mesh,toto,nodeindexes,matid,index);
// }
//template<> class TPZGeoElRefPattern<pzgeom::TPZBlendNACA>;


template class
TPZRestoreClass< TPZGeoElRefPattern<TPZBlendNACA>>;


