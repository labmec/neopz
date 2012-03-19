// $Id: TPZYCVonMises.h,v 1.12 2009-06-29 22:54:01 erick Exp $

#ifndef TPZYCVONMISES_H
#define TPZYCVONMISES_H

#include "TPZTensor.h"
#include "pzfmatrix.h"

#include "pzlog.h"

/**
Implementa  a plastificacao do criterio de Von Mises
*/
class TPZYCVonMises {
    

public:

  enum {NYield = 1};
	
    const char * Name() const
    {
	   return "TPZYCVonMises";	
    }
	
    void Print(std::ostream & out) const
    {
       out << Name();
    }
	
	int GetForceYield()
	{
		return 0; // nothing to be done in this yield criterium
	}
	
	void SetForceYield(const int forceYield)
	{
		// nothing to be done in this yield criterium
	}

	/**
	 * Checks if the proposed yield state leads to post-peak material behaviour. If so, the material
	 * is forced to behave in post-peak in order to avoid equation switching during Newton's method
	 * in the PlasticLoop routines.
	 * @param [in] sigma stress state
	 * @param [in] A Thermo Force
	 */
	void SetYieldStatusMode(const TPZTensor<REAL> & sigma, const REAL & A)
	{
		// nothing to be done in this yield criterium
	}
	
    /**
    Calculo do criterio de plastificacao 
    @param [in] sigma tensao atual
    @param [in] A forca thermodinamica atual
	@param [in] checkForcedYield indicates wether to force post-peak failure behavior
    */  
    template < class T>
    void Compute(const TPZTensor<T> & sigma, const T & A, TPZVec<T> &res, int checkForcedYield = 0) const;
    
    /**
    Derivada da funcao de plastificacao
    @param [in] sigma tensao atual
    @param [in] A forca termodinamica atual
    @param [out] Derivida com respeito a tensao
	@param [in] checkForcedYield indicates wether to force post-peak failure behavior
    */
    template <class T> 
    void N(const TPZTensor<T> & sigma,const T & A,  TPZVec<TPZTensor<T> > & Ndir, int checkForcedYield = 0) const;

    /**
    Derivada da funcao de plastificacao com respeito a forca termodinamica
    @param [in] sigma tensao atual
    @param [in] A forca termodinamica atual
    @param [out] Derivida com respeito a forca termodinamica
    @param [in] checkForcedYield indicates wether to force post-peak failure behavior
    */
    template <class T> 
    void H(const TPZTensor<T> & sigma,const T & A,  TPZVec<T> & h, int checkForcedYield = 0) const;

public:
//////////////////CheckConv related methods/////////////////////

    /**
    number of types of residuals
    */
    int NumCases()
    {
      return 3;
    }

    static TPZTensor<REAL> gRefTension;
    /**
    LoadState will keep a given state as static variable of the class
    */
    void LoadState(TPZFMatrix<REAL> &state)
    {
    #ifdef LOG4CXX_PLASTICITY
        LoggerPtr logger(Logger::getLogger("plasticity.ycvonmises"));
    #endif
      int i;
      for(i=0; i<6; i++) gRefTension.fData[i] = state(i,0);
	#ifdef LOG4CXX_PLASTICITY
      std::stringstream sout;
      sout << "Tensao " << state;
      LOGPZ_DEBUG(logger,sout.str().c_str());
	#endif
    }

    void ComputeTangent(TPZFMatrix<REAL> &tangent, TPZVec<REAL> &, int icase)
    {
    #ifdef LOG4CXX_PLASTICITY
        LoggerPtr logger(Logger::getLogger("plasticity.ycvonmises"));
    #endif
      TPZTensor<REAL> dj2;
      REAL A(1.e7);
      TPZVec<REAL> phivec(1,0.);
      TPZVec<TPZTensor<REAL> > ndir(1);
      int i;
      REAL j2, mult;
      switch(icase)
      {
        case 0:
          gRefTension.dJ2(dj2);
          tangent.Redim(1,6);
          for(i=0; i<6; i++) tangent(0,i) = dj2.fData[i];
          break;
        case 1:
          j2 = gRefTension.J2();
          gRefTension.dJ2(dj2);
          tangent.Redim(1,6);
          mult = 0.5/sqrt(j2);
          dj2.Multiply(mult,1.);
          for(i=0; i<6; i++) tangent(0,i) = dj2.fData[i];
          break;
        case 2:
          N(gRefTension,A,ndir,0);
          tangent.Redim(1,6);
          for(i=0; i<6; i++) tangent(0,i) = ndir[0].fData[i];
      }
	#ifdef LOG4CXX_PLASTICITY
      std::stringstream sout;
      sout << "Matriz tangent " << tangent;
      LOGPZ_DEBUG(logger,sout.str().c_str());
	#endif
    }

    void Residual(TPZFMatrix<REAL> &res,int icase)
    {
    #ifdef LOG4CXX_PLASTICITY
        LoggerPtr logger(Logger::getLogger("plasticity.ycvonmises"));
    #endif
      TPZTensor<REAL> gradtheta;
      REAL A(1.e7),j2;
      TPZVec<REAL> phivec(1,0.);
      switch(icase)
      {
        case 0:
          j2 = gRefTension.J2();
          res.Redim(1,1);
          res(0,0) = j2;
          break;
        case 1:
          j2 = gRefTension.J2();
          res.Redim(1,1);
          res(0,0) = sqrt(j2);
          break;
        case 2:
          Compute(gRefTension,A,phivec,0);
          res.Redim(1,1);
          res(0,0) = phivec[0];
          break;
      }
	#ifdef LOG4CXX_PLASTICITY
      std::stringstream sout;
      sout << "valor phi " << res(0,0);
      LOGPZ_DEBUG(logger,sout.str().c_str());
	#endif
    }

//////////////////CheckConv related methods/////////////////////

};

template <class T>
void TPZYCVonMises::Compute(const TPZTensor<T> & sigma,const T & A, TPZVec<T> &res, int checkForcedYield) const{
	
	T J2 = sigma.J2();
	T temp = sqrt(T(3.)*J2);
	T result = temp-A;
    res[0] = result;
}

template <class T> 
void TPZYCVonMises::N(const TPZTensor<T> & sigma,const T & A, TPZVec<TPZTensor<T> > & Ndir, int checkForcedYield) const
{
	T  temp1= sqrt(T(3.)/T(2.));
	TPZTensor<T> s;
	sigma.S(s);
	T NORM	= sqrt(s.XX()*s.XX()+s.YY()*s.YY()+s.ZZ()*s.ZZ()+T(2.)*(s.XY()*s.XY())+T(2.)*(s.XZ()*s.XZ())+T(2.)*(s.YZ()*s.YZ()));
	s.Multiply(T(1.)/NORM,T(1.));
	s.Multiply(temp1, T(1.));
	Ndir[0]=s;
//	sigma.dJ2(Ndir[0]);
//	Ndir[0].Multiply(T(0.5)/sqrt(sigma.J2()),sqrt(3.));
	
}

template <class T> 
void TPZYCVonMises::H(const TPZTensor<T> & sigma,const T & A, TPZVec<T> & h, int checkForcedYield) const
{
    h[0] = 1.;
// Conflicting definitions:
// Eduardo Souza Neto, pg 242 -> H = -1.0 for von Mises
//                               H =  dKsi/dA
//                     pg 146 -> H = -dKsi/dA

}

#endif //TPZYCVONMISES_H
