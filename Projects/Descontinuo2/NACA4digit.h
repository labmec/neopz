
#ifndef TPZNACA4DIGITSH
#define TPZNACA4DIGITSH

#include <math.h>
#include "pzreal.h"
#include "pzvec.h"

class TPZNACAXXXX {
  double fCord;
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

  double P()
  {
    int aux = fFourDigits/100;
    aux -= ((int)(aux/10))*10;
    return (double)aux/10.;
}

double M()
{
   int aux = fFourDigits/1000;
   return (double)aux/100.*fCord;
}

double TT()
{
   int aux = fFourDigits - ((int)(fFourDigits/100))*100;
   return (double)aux/100.*fCord;
}

// Mean line for the wing
double yc(double x)
{
   if(x/fCord<fP)
   {
      return fM/fP/fP*(2.*fP*x/fCord-x*x/fCord/fCord);
   }else
   {
      return fM/(1.-fP)/(1.-fP)*(1.-2.*fP+2.*fP*x/fCord-x*x/fCord/fCord);
   }
}

double dyc(double x)
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
double yt(double x)
{
   double aux = x/fCord;
   const double a0 = 1.4845,
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
double xu(double x)
{
   return x-yt(x)*sin(atan(dyc(x)));
}

double yu(double x)
{
   return yc(x) + yt(x)*cos(atan(dyc(x)));
}


// inferior profile
double xl(double x)
{
   return x+yt(x)*sin(atan(dyc(x)));
}

double yl(double x)
{
   return yc(x) - yt(x)*cos(atan(dyc(x)));
}

  REAL NearestParameter(TPZVec<REAL> &pt, int &upperlower);

 public:

// with attack angle
// superior profile
double xua(double x)
{
   return fX0[0]+(xu(x)-fCord/2.)*cos(fAngle) + yu(x) * sin(fAngle) + fCord/2.;
}

double yua(double x)
{
   return fX0[1]+yu(x)*cos(fAngle) - (xu(x)-fCord/2.) * sin(fAngle);
}

// inferior profile
double xla(double x)
{
   return fX0[0]+(xl(x)-fCord/2.)*cos(fAngle) + yl(x) * sin(fAngle) + fCord/2.;
}

double yla(double x)
{
   return fX0[1]+yl(x)*cos(fAngle) - (xl(x)-fCord/2.) * sin(fAngle);
}

  void ProjectPoint(TPZVec<REAL> &pt) 
  {
    int uplow;
    REAL par = NearestParameter(pt,uplow);
    if(uplow == 0) {
      pt[0] = xla(par);
      pt[1] = yla(par);
    } else {
      pt[0] = xua(par);
      pt[1] = yua(par);
    }
  }

};

#endif //NACA4DIGITSH
