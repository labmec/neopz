// $Id: TPZYCRankine.h,v 1.2 2010-06-11 22:12:14 diogo Exp $

#ifndef TPZYCRANKINE_H
#define TPZYCRANKINE_H

#include "TPZTensor.h"
#include "pzlog.h"

template <class T_YCBASE>
class TPZYCRankine : public T_YCBASE{

public:
  
    TPZYCRankine():T_YCBASE(),fYieldT(0.)
	{
	}
	
	~TPZYCRankine()
	{
	}
	
	enum {NYield = 1 + T_YCBASE::NYield};
    
    const char * Name() const
    {
	   return "TPZYCRankine<>";	
    }
	    		
	/**
	 * Sets the material tensile strength
	 * @param [in] YieldT tensile strength
	 */
	void SetUpRankine(const REAL & YieldT)
	{
		fYieldT = YieldT;
	}
	
  /**
   * Evaluate the yield criteria
   * @param [in] sigma current stress tensor
   * @param [in] A current thermodynamical force
   * @param [out] result
   * @param [in] checkForcedYield indicates wether to force post-peak failure behavior
   */
  template < class T>
  void Compute(const TPZTensor<T> & sigma, const T & A, TPZVec<T> &result, int checkForcedYield = 0) const;

  /**
   * Derivative of the yield function
   * @param [in] sigma current stress tensor
   * @param [in] A current thermodynamical force
   * @param [out] Ndir Stress derivative
   * @param [in] checkForcedYield indicates wether to force post-peak failure behavior
   */
  template <class T>
  void N(const TPZTensor<T> & sigma,const T & A,  TPZVec<TPZTensor<T> > & Ndir, int checkForcedYield = 0) const;

  /**
   * Derivative of the yield function with respect to the thermodynamical force
   * @param [in] sigma current stress tensor
   * @param [in] A current thermodynamical force
   * @param [out] h Derivative with respect to thermodynamical force
   * @param [in] checkForcedYield indicates wether to force post-peak failure behavior
   */
  template <class T>
	void H(const TPZTensor<T> & sigma,const T & A,  TPZVec<T> & h, int checkForcedYield = 0) const;

protected:
	
	REAL fYieldT;

};

/**
 * Calculo do criterio de plastificacao
 * @param [in] sigma tensao atual
 * @param [in] A forca thermodinamica atual
 */
template < class T_YCBASE>
template < class T>
inline void TPZYCRankine<T_YCBASE>::Compute(const TPZTensor<T> & sigma, const T & A,TPZVec<T> &result, int checkForcedYield) const
{
	T_YCBASE::Compute(sigma, A, result, checkForcedYield);
	// Chen Plasticity for structural engineers pg 88
	T theta, I1, J2, J3;
	I1 = sigma.I1();
	J2 = sigma.J2();
	J3 = sigma.J3();
	
	T sqrtJ2 = sqrt(J2);
	T sqrtJ232 = J2 * sqrtJ2;
	//T value = J3 / pow(J2, 3./2.) * T(sqrt(27.)/2.);
	//T value = J3 / exp ( T(3./2.) * log(J2) ) * T(sqrt(27.)/2.);
	T value = J3 / sqrtJ232 * T(sqrt(27.)/2.);
	theta = acos( value / T(3.) );
	// theta MUST lie between 0 and 60 deg !!!
	
	value = sqrt(J2) * T(2. * sqrt(3.)) * cos(theta) + I1 - T(3. * fYieldT);
	result[T_YCBASE::NYield] = value;
}


/**
 * Derivada da funcao de plastificacao
 * @param [in] sigma tensao atual
 * @param [in] A forca termodinamica atual
 * @param [out] Derivida com respeito a tensao
 */
template < class T_YCBASE>
template <class T> 
inline void TPZYCRankine<T_YCBASE>::N(const TPZTensor<T> & sigma,const T & A,  TPZVec<TPZTensor<T> > & Ndir, int checkForcedYield) const
{
	T_YCBASE::N(sigma, A, Ndir, checkForcedYield);
	// implementar contribuicao para 
	
	TPZTensor<T> NRankine;
	NRankine.XX() = T(0.);//dfdsigmaxx
	NRankine.YY() = T(0.);//dfdsigmayy
	NRankine.ZZ() = T(0.);//dfdsigmazz
	NRankine.XY() = T(0.);//dfdsigmaxy
	NRankine.XZ() = T(0.);//dfdsigmaxz
	NRankine.YZ() = T(0.);//dfdsigmayz
	
	Ndir[T_YCBASE::NYield] = NRankine;
}

template < class T_YCBASE>
template <class T> 
inline void TPZYCRankine<T_YCBASE>:: H(const TPZTensor<T> & sigma,const T & A,  TPZVec<T> & h, int checkForcedYield)const
{
	h[T_YCBASE::NYield] = 0.; // not really necessary since this model does not handle hardening.
}

#endif //TPZYCRANKINE_H
