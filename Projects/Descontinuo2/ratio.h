#ifndef RATIO_H
#define RATIO_H

// used by xpg
double Spg(double a0, double q, int N)
{
   return a0*(pow(q, N+1.)-1.)/(q-1.);
}

// returns a real between 0 and 1, as the
// percentage of the length according to the
// q parameter. (ratio)
double xpg(double q, int n, int N)
{
   if(q==1.)return ((double)n) / ((double)N);
   return Spg(1., q, n-1)/Spg(1., q, N-1);
}

double PI = 3.14159265359;

// return a real between 0 and 1, but with
// a trigonometric dispersion.
// Values are concentrated in the extremes.
double xtrig(int n, int N)
{
   return (1.-cos(((double)n)/((double)N)*PI))/2.;
}

#endif
