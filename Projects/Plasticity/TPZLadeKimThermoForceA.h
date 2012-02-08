// $Id: TPZLadeKimThermoForceA.h,v 1.10 2009-06-08 03:19:34 erick Exp $

#ifndef TPZLADEKIMTHERMOFORCEA_H
#define TPZLADEKIMTHERMOFORCEA_H

#include "pzvec.h"
#include "pzreal.h"

#ifndef CHECKCONV
#define CHECKCONV
#include "checkconv.h"
#endif
/**
Classe que implementa o calculo da forca termodinamica (Souza Neto p. 144) e suas derivadas
Neste caso utiliza-se encruamento linear
*/
class TPZLadeKimThermoForceA {
public:

    TPZLadeKimThermoForceA() : fRho(0), fD(0), fPa(0)
    {
    }
	
    TPZLadeKimThermoForceA(const TPZLadeKimThermoForceA & source)
    {
	   fRho	= source.fRho;
       fD	= source.fD;
	   fPa	= source.fPa;
    }

    TPZLadeKimThermoForceA & operator=(const TPZLadeKimThermoForceA & source)
    {
	   fRho	= source.fRho;
       fD	= source.fD;
	   fPa	= source.fPa;
		return *this;
    }

	const char * Name() const
    {
	   return "TPZLadeKimThermoForceA";	
    }
	
    void Print(std::ostream & out) const
    {
		out << "\n" << this->Name();
		out << "\n fRho  = " << fRho;
		out << "\n fD     = " << fD;
		out << "\n fPa    = " << fPa;
    }


	
    void SetUp(REAL ksi1, REAL p, REAL h, REAL C, REAL pa) 
    {
        fRho = p / h;
        fD   = C / pow( 27. * ksi1 + 3., fRho);
        fPa  = pa;
    }


    /**
    Calculo do valor da forca termo dinamica
    */
    template <class T>
    T Compute(const T & alpha) const;

    /**
    * Specialization
    */
    REAL Compute(const REAL & alpha) const;

    /**
    Calculo da derivada da forca termodinamica
    */
    template <class T>
    T ComputeTangent(const T & alpha) const;

public:

   /**
     Parameter related to the ThermoForceA (Hardening)
     Rho is the inverse of the exponent of the plastic work in the hardening function.
     It is compute from p and h.
     ->p defines the exponent of the evolution of plastic work along with I1 in pure isotropic
     compression. (Wp = C * I1^p)
     ->the h constant (YC) models the curvature of the YC meridians, i.e., how the meridians
     vary along with the level of hydrostatic stress (I1/3).
   */
   REAL fRho;

   /**
     Parameter related to the ThermoForceA (Hardening)
     D is the constant multiplier in the hardening function.
     It is computed from Ksi1, C and Rho
     ->C defines the constant of the evolution of plastic work along with I1 in pure isotropic
     compression. (Wp = C * I1^p)
     ->Ksi1 models the shape of the yield funcition at the deviatoric planes.
     ->rho is defined above.
   */
   REAL fD;

   /**
     Atmospheric pression to hold stress conversion constant
   */
   REAL fPa;


public:
//////////////////CheckConv related methods/////////////////////

    /**
    number of types of residuals
    */
    int NumCases()
    {
      return 1;
    }

    static REAL gRefThermoForce;
    /**
    LoadState will keep a given state as static variable of the class
    */
    void LoadState(TPZFMatrix &state)
    {
    #ifdef LOG4CXX_PLASTICITY
        LoggerPtr logger(Logger::getLogger("plasticity.ladekimthermoforce"));
    #endif
      gRefThermoForce = state(0,0);
	#ifdef LOG4CXX_PLASTICITY
      std::stringstream sout;
      sout << "Plastic Work " << state;
      LOGPZ_DEBUG(logger,sout.str().c_str());
	#endif
    }

    void ComputeTangent(TPZFMatrix &tangent, TPZVec<REAL> &, int icase)
    {
    #ifdef LOG4CXX_PLASTICITY
        LoggerPtr logger(Logger::getLogger("plasticity.ladekimthermoforce"));
    #endif

      switch(icase)
      {
        case 0:
          //Compute
          tangent.Redim(1,1);
          tangent(0,0) = ComputeTangent(gRefThermoForce);
          break;

          break;
      }
	#ifdef LOG4CXX_PLASTICITY
      std::stringstream sout;
      sout << "Matriz tangent " << tangent;
      LOGPZ_DEBUG(logger,sout.str().c_str());
    #endif
    }

    void Residual(TPZFMatrix &res,int icase)
    {
    #ifdef LOG4CXX_PLASTICITY
        LoggerPtr logger(Logger::getLogger("plasticity.ladekimthermoforce"));
    #endif
      switch(icase)
      {
        case 0:
          //Compute
          res.Redim(1,1);
          res(0,0) = Compute(gRefThermoForce);
        break;

      }
    #ifdef LOG4CXX_PLASTICITY
      std::stringstream sout;
      sout << "Residual vector " << res;
      LOGPZ_DEBUG(logger,sout.str().c_str());
    #endif
    }

static void CheckConv()
{
   // Creating the LadeNelsonElasticResponse obejct
   TPZLadeKimThermoForceA TFALadeNelson;
   // setup with data from Sacramento River Sand
   // Lade, Poul V., Kim, Moon K. Single Hardening Constitutive Model
   // for soil, Rock and Concrete. Int. J. Solids Structures, vol32 
   // No14 pp 1963-1978,1995
   REAL m = 0.093,
        p = 1.82,
        h = 0.765,
        C = 0.396e-4,
        pa = 14.7;
   REAL ksi1 = 0.00155 * pow(m, -1.27);
   TFALadeNelson.SetUp(ksi1, p, h, C, pa);
   TPZFNMatrix<1> Alpha(1, 1, 1.7);
   TPZFNMatrix<1> Range(1, 1, Alpha(0,0)*(1./19.));
   TPZVec< REAL > Coefs(1,1.);
   CheckConvergence(TFALadeNelson, Alpha, Range, Coefs);
}

//////////////////CheckConv related methods/////////////////////

};


template < class T >
inline T TPZLadeKimThermoForceA::Compute(const T & alpha) const
{
    T localAlpha(alpha);
    //if(fabs(shapeFAD::val(localAlpha) ) < 1.e-60)localAlpha.val()+=1.e-60;
	if(shapeFAD::val(localAlpha) < 1.e-60)localAlpha.val() = 1.e-60;
    REAL invRho = 1./fRho;
    //return T( pow( fD, - invRho ) ) * pow( localAlpha / T(fPa), invRho );
	return T( pow( fD, - invRho ) ) * exp( T(invRho) * log( localAlpha / T(fPa) ) );
}

inline REAL TPZLadeKimThermoForceA::Compute(const REAL & alpha) const
{
    REAL localAlpha = alpha;
    //if(fabs(localAlpha) < 1.e-60)localAlpha+=1.e-60;
	if( localAlpha < 1.e-60)localAlpha = 1.e-60;
    //if(fabs(alpha) < 1.e-12) return 0.;
    REAL invRho = 1./fRho;
    return pow( fD, - invRho ) * pow( alpha / fPa, invRho );
}

template < class T >
inline T TPZLadeKimThermoForceA::ComputeTangent(const T & alpha) const
{
    REAL invRho = 1./fRho;
    return T( pow( fD, - invRho ) * invRho / fPa ) * pow( alpha / fPa, invRho-1.);
}

#endif //TPZTHERMOFORCEA_H
