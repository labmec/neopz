// $Id: TPZElasticResponse.h,v 1.8 2009-07-04 03:31:24 erick Exp $

#ifndef TPZELASTICRESPONSE_H
#define TPZELASTICRESPONSE_H

#include "TPZTensor.h"
#include "pzreal.h"

/**
Calcula a tensao em funcao de deformacao (elastica)
*/
class TPZElasticResponse {
public:

   TPZElasticResponse() : fLambda(0.), fMu(0.)
   {
   }
	
   TPZElasticResponse(const TPZElasticResponse & source)
   {
	   fLambda	= source.fLambda;
       fMu		= source.fMu;
   }

   TPZElasticResponse & operator=(const TPZElasticResponse & source)
   {
	   fLambda	= source.fLambda;
       fMu		= source.fMu;
	   return *this;
   }
	
   /**
   Construtor da classe em funcao das variaveis de lame
   @param [in] lambda primeira variavel de lame
   @param [in] mu segunda variavel de lame
   */
   void SetUp(REAL elast, REAL poisson)
   {
   	
   	fLambda = poisson*elast/((1.+poisson)*(1.-2.*poisson));
	fMu = elast/(2.*(1.+poisson));
   }
	
   const char * Name() const
   {
	   return "TPZElasticResponse";	
   }
	
   void Print(std::ostream & out) const
   {
       out << this->Name();
	   out << "\n fLambda = " << fLambda;
	   out << "\n fMu = " << fMu;
   }
	
   /**
   Calcula o tensor de tensao em funcao do tensor de deformacao
   @param [in] epsilon tensor de deformacao
   @param [out] sigma tensor de tensao
   */
    template <class T>
    void Compute(TPZTensor<T> & epsilon, TPZTensor<T> & sigma) const
    {
    T trace = epsilon.I1();
    sigma.Identity();
    sigma.Multiply(trace,fLambda);
    sigma.Add(epsilon,2.*fMu);
    }
	
	/**
	 Computes the elastic matrix, writing it to the matrix Kep 
	 @param Kef[out] matrix to write to
	*/
    void ElasticMat(TPZFMatrix & Kef)
	{
		REAL Mu2    = 2 * fMu;
		
		Kef.Zero();
		
		Kef(_XX_, _XX_) += fLambda;
		Kef(_XX_, _YY_) += fLambda;
		Kef(_XX_, _ZZ_) += fLambda;
		Kef(_YY_, _XX_) += fLambda;
		Kef(_YY_, _YY_) += fLambda;
		Kef(_YY_, _ZZ_) += fLambda;
		Kef(_ZZ_, _XX_) += fLambda;
		Kef(_ZZ_, _YY_) += fLambda;
		Kef(_ZZ_, _ZZ_) += fLambda;
		
		int i;
		for(i = 0; i < 6; i++)Kef(i, i) += Mu2;
	}
    
    REAL fLambda;
    REAL fMu;
};


#endif //TPZELASTICRESPONSE_H

