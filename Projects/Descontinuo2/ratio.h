#ifndef RATIO_H
#define RATIO_H

#include "pzreal.h"

#include <cmath>

using namespace std;

// used by xpg
inline REAL Spg(REAL a0, REAL q, int N)
{
   return a0*(pow(q, (REAL)(N+1.))-1.)/(q-1.);
}

// returns a real between 0 and 1, as the
// percentage of the length according to the
// q parameter. (ratio)
inline REAL xpg(REAL q, int n, int N)
{
   if(q==1.)return ((REAL)n) / ((REAL)N);
   return Spg(1., q, n-1)/Spg(1., q, N-1);
}


// return a real between 0 and 1, but with
// a trigonometric dispersion.
// Values are concentrated in the extremes.
inline REAL xtrig(int n, int N)
{
   return (1.-cos(((REAL)n)/((REAL)N)*M_PI))/2.;
}

#endif
