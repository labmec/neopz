#ifndef FADTYPE_H
#define FADTYPE_H



#include "fad.h"
#include "tinyfad.h"
#include "tfad.h" // tinyFadET


#include "pzreal.h"
#include "pzvec.h"

const int REALdim = 3;

typedef TFad<REALdim, REAL> TFADETREAL_;//  TinyFAD with Expression templates
typedef TinyFad<REALdim, REAL> TFADREAL_;// TinyFAD - the easiest to use
typedef Fad<REAL> FADREAL_; // the simplest FAD

typedef FADREAL_ FADREAL;

typedef Fad< Fad<REAL> > FADFADREAL;
typedef Fad< Fad<STATE> > FADFADSTATE;


#include "pzfmatrix.h"


class shapeFAD
{
public:
    static void ExplodeDerivatives(TPZVec<FADREAL> & in, TPZFMatrix<REAL> & phiOut, TPZFMatrix<REAL> & dphiOut)
    {
        int nphi = in.NElements(), nder = in[0].size();
        phiOut.Redim(nphi, 1);
        dphiOut.Redim(nder, nphi);
	
        int i, j;
		for(i = 0; i < nphi; i++)
        {
            phiOut(i,0) = in[i].val();
            for(j = 0; j < nder; j++)
            {
                dphiOut(j,i) = in[i].dx(j);
            }
        }
	    return;
    }

    static REAL val( const int number)
    {
        return (REAL)number;
    }
    static REAL val( const int64_t number)
    {
        return (REAL)number;
    }
    static REAL val( const float number) 
    {
        return (REAL)number;
    }
    static REAL val( const double number) 
    {
        return (REAL)number;
    }
    static REAL val( const long double number) 
    {
        return (REAL)number;
    }
    static REAL val( const std::complex<float> number) 
    {
        return (REAL)number.real();
    }
    static REAL val( const std::complex<double> number) 
    {
        return (REAL)number.real();
    }
    static REAL val( const std::complex<long double> number) 
    {
        return (REAL)number.real();
    }

    template <class T>
    static REAL val(const T number) 
    {
        return TPZExtractVal::val( number.val() ); // recursively downgrading until REAL type is reached
    }

    template<class T>
    static bool IsZero( const T & a ){
      return ::IsZero(TPZExtractVal::val(a));
    }

};

template class Fad<float>;
template class Fad<double>;
template class Fad<long double>;


#endif

