
#ifndef TPZNACA4DIGITSH
#define TPZNACA4DIGITSH

#include <math.h>
#include "pzreal.h"
#include "pzvec.h"

class TPZNACAXXXX {
  REAL fCord;
  int fFourDigits;
  REAL fAngle;
  REAL fX0[3];
  REAL fP;
  REAL fM;
  REAL fTT;
  
public:
  
  TPZNACAXXXX(REAL cord, int FourDigits, REAL angle, TPZVec<REAL> &x0) :
    fCord(cord), fFourDigits(FourDigits), fAngle(angle)
  {
    fX0[0] = x0[0];
    fX0[1] = x0[1];
    fX0[2] = x0[2];
    fP = P();
    fM = M();
    fTT = TT();
  }
  

 private:

  REAL P()
  {
    int aux = fFourDigits/100;
    aux -= ((int)(aux/10))*10;
    return (REAL)(aux/10.);
}

REAL M()
{
   int aux = fFourDigits/1000;
   return (REAL)(aux/100.)*fCord;
}

REAL TT()
{
   int aux = fFourDigits - ((int)(fFourDigits/100))*100;
   return (REAL)(aux/100.)*fCord;
}

// Mean line for the wing
REAL yc(REAL x)
{
   if(x/fCord<fP)
   {
      return fM/fP/fP*(2.*fP*x/fCord-x*x/fCord/fCord);
   }else
   {
      return fM/(1.-fP)/(1.-fP)*(1.-2.*fP+2.*fP*x/fCord-x*x/fCord/fCord);
   }
}

REAL dyc(REAL x)
{
   if(x/fCord<fP)
   {
      return 2.*fM/fP/fP*(fP-x/fCord)/fCord;
   }else
   {
      return 2.*fM/(1.-fP)/(1.-fP)*(fP-x/fCord)/fCord;
   }
}

// thickness
REAL yt(REAL x)
{
   REAL aux = x/fCord;
   const REAL a0 = 1.4845,
                a1 = -.6300,
		a2 = -1.7580,
		a3 = 1.4215,
		a4 = -.5075;

   return fTT * (a0*sqrt(aux)+
                 a1*aux+
		 a2*aux*aux+
		 a3*aux*aux*aux+
		 a4*aux*aux*aux*aux);
}

// superior profile
REAL xu(REAL x)
{
   return x-yt(x)*sin(atan(dyc(x)));
}

REAL yu(REAL x)
{
   return yc(x) + yt(x)*cos(atan(dyc(x)));
}


// inferior profile
REAL xl(REAL x)
{
   return x+yt(x)*sin(atan(dyc(x)));
}

REAL yl(REAL x)
{
   return yc(x) - yt(x)*cos(atan(dyc(x)));
}

  REAL NearestParameter(TPZVec<REAL> &pt, int &upperlower, int maxPt);

 public:

// with attack angle
// superior profile
REAL xua(REAL x)
{
   return fX0[0]+(xu(x)-fCord/2.)*cos(fAngle) + yu(x) * sin(fAngle) + fCord/2.;
}

REAL yua(REAL x)
{
   return fX0[1]+yu(x)*cos(fAngle) - (xu(x)-fCord/2.) * sin(fAngle);
}

// inferior profile
REAL xla(REAL x)
{
   return fX0[0]+(xl(x)-fCord/2.)*cos(fAngle) + yl(x) * sin(fAngle) + fCord/2.;
}

REAL yla(REAL x)
{
   return fX0[1]+yl(x)*cos(fAngle) - (xl(x)-fCord/2.) * sin(fAngle);
}

  void ProjectPoint(TPZVec<REAL> &pt, int maxPt = 1000)
  {
    int uplow;
    REAL par = NearestParameter(pt,uplow, maxPt);
    if(uplow == 0) {
      pt[0] = xla(par);
      pt[1] = yla(par);
    } else {
      pt[0] = xua(par);
      pt[1] = yua(par);
    }
  }

};

#include "pzreal.h"
#include "pzvec.h"
#include "pzeltype.h"
#include "pzgmesh.h"
#include "pzflowcmesh.h"
#include "pzeulerconslaw.h"

void NACAPoints(TPZNACAXXXX &profile, TPZVec< TPZVec<REAL> > & pt, TPZVec< TPZVec< int> > &elms, int nSubdiv);
TPZGeoMesh * CreateNACAGeoMesh(TPZGeoMesh *gmesh, TPZNACAXXXX &profile, TPZVec< TPZVec< REAL > > & nodes,
                               TPZVec< TPZVec< int > > & elms,
                               MElementType ElType, int matId,
                               TPZVec<TPZGeoEl *> & gEls,
                               int nSubdiv);
// Creating all the geometric and computational meshes
// for the reflected shock problem.

TPZFlowCompMesh *
NACACompMesh(TPZFlowCompMesh *cmesh, REAL CFL, REAL delta,
             int degree, int nSubdiv,
             TPZArtDiffType DiffType,
             TPZTimeDiscr Diff_TD,
             TPZTimeDiscr ConvVol_TD,
             TPZTimeDiscr ConvFace_TD, std::ostream &options);


#endif //NACA4DIGITSH
