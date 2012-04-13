/**
 * Numerical Recipes - The art of Scientific Computing
 * 3rd Edition
 *
 * William H. Press
 * Saul A. Teukolsky
 * William T. Vetterling
 * Brian P. Flannery
 */

#ifndef ADAPTH
#define ADAPTH

#include "pzreal.h"
#include "limits"

struct Adapt
{

public:
	
	Adapt(REAL tol);
    
	template <class T>
	REAL integrate(T &func, const REAL a, const REAL b);



private:
    
    template <class T>
	REAL adaptlob(T &func, const REAL a, const REAL b, const REAL fa, const REAL fb, const REAL is);
    
    REAL TOL,toler;
    bool terminate,out_of_tolerance;
    static const REAL alpha,beta,x1,x2,x3,x[12];
};


inline Adapt::Adapt(REAL tol) : TOL(tol),terminate(true),out_of_tolerance(false)
{
	const REAL EPS = std::numeric_limits<REAL>::epsilon();
	if(TOL < 10.0 * EPS)
    {
		TOL = 10.0 * EPS;
    }
}

template <class T>
inline REAL Adapt::integrate(T &func, const REAL a, const REAL b)
{
	REAL m,h,fa,fb,i1,i2,is,erri1,erri2,r,y[13];
	m = 0.5 * (a + b);
	h = 0.5 * (b-a);
	fa = y[0] = func(a);
	fb = y[12] = func(b);
	for(int i = 1; i < 12; i++ )
    {
		y[i] = func(m + x[i] * h);
    }
	i2 = (h/6.0) * (y[0] + y[12] + 5.0 * (y[4] + y[8]));
	i1 = (h/1470.0) * (77.0 * (y[0] + y[12]) + 432.0 * (y[2] + y[10]) + 625.0 * (y[4] + y[8]) + 672.0 * y[6]);
	is = h * (0.0158271919734802 * (y[0] + y[12]) + 0.0942738402188500 * (y[1] + y[11]) + 0.155071987336585 * (y[2] + y[10])  + 
              0.188821573960182 * (y[3] + y[9]) + 0.199773405226859 * (y[4] + y[8]) + 0.224926465333340 * (y[5] + y[7]) + 0.242611071901408 * y[6]);
	erri1 = fabs(i1 - is);
	erri2 = fabs(i2 - is);
	r = (erri2 !=  0.0) ? erri1/erri2 : 1.0;
	toler = (r > 0.0 && r < 1.0) ? TOL/r : TOL;
	if(is == 0.0)
    {
		is = b-a;
    }
	is = fabs(is);
    
    REAL answ = adaptlob(func,a,b,fa,fb,is);
    
	return answ;
}

template <class T>
inline REAL Adapt::adaptlob(T &func, const REAL a, const REAL b, const REAL fa, const REAL fb, const REAL is)
{
	REAL m,h,mll,ml,mr,mrr,fmll,fml,fm,fmrr,fmr,i1,i2;
	m = 0.5 * (a + b);
	h = 0.5 * (b-a);
	mll = m-alpha * h;
	ml = m-beta * h;
	mr = m + beta * h;
	mrr = m + alpha * h;
	fmll = func(mll);
	fml = func(ml);
	fm = func(m);
	fmr = func(mr);
	fmrr = func(mrr);
	i2 = h/6.0 * (fa + fb + 5.0 * (fml + fmr));
	i1 = h/1470.0 * (77.0 * (fa + fb) + 432.0 * (fmll + fmrr) + 625.0 * (fml + fmr) + 672.0 * fm);
	if(fabs(i1 - i2) <=  toler * is || mll <=  a || b <=  mrr) 
    {
		if((mll <=  a || b <=  mrr) && terminate) 
        {
			out_of_tolerance = true;
			terminate = false;
		}
		return i1;
	}
	else
    {
        REAL val =  adaptlob(func,a,mll,fa,fmll,is) + adaptlob(func,mll,ml,fmll,fml,is) + adaptlob(func,ml,m,fml,fm,is) + 
        adaptlob(func,m,mr,fm,fmr,is) + adaptlob(func,mr,mrr,fmr,fmrr,is) + adaptlob(func,mrr,b,fmrr,fb,is);
        
		return val;
    }
}


#endif
