/**
 * @file
 */

#ifndef TPZYCTRESCAREGULARIZED_H
#define TPZYCTRESCAREGULARIZED_H

#include "TPZYCTresca.h"
#include "pzlog.h"
#include <math.h>
#include <iostream> 

/**
 * @brief This class implements a tresca yield criteria whrere the gradient of the inverse angle is regularized.
 * @author LabMeC
 */
class TPZYCTrescaRegularized : public TPZYCTresca
{
public:
	
    public:
int ClassId() const override;

    
    const char * Name() const
    {
		return "TPZYCTrescaRegularized";	
    }
	
    void Print(std::ostream & out) const override
    {
		out << Name();
    }
	
	/**
	 * @brief Calculo do criterio de plastificacao
	 * @param[in] sigma tensao atual
	 * @param[in] A forca thermodinamica atual
	 * @param[out] result Derivative
	 * @param[in] checkForcedYield indicates wether to force post-peak failure behavior
	 */	
	template < class T>
	void Compute(const TPZTensor<T> & sigma, const T & A,TPZVec<T> &result, int checkForcedYield = 0) const;
	
protected:
	/**
	 * @brief Compute the inverse angle of the tresca yield criterium formula and
	 * the related data
	 * @param[in] sigma stress tensor
	 * @param[out] theta one third of the asin of the inverse angle
	 * @param[out] gradasin gradient of the inverse angle
	 */
	template <class T>
	void GradTheta(const TPZTensor<T> & sigma,T & theta, TPZTensor<T> & gradasin) const;
	
public:
	/**
	 * @brief Derivada da funcao de plastificacao
	 * @param[in] sigma tensao atual
	 * @param[in] A forca termodinamica atual
	 * @param[out] Ndir Derivada com respeito a tensao
	 * @param[in] checkForcedYield indicates wether to force post-peak failure behavior
	 */
	template <class T> 
	void N(const TPZTensor<T> & sigma,const T & A,  TPZVec<TPZTensor<T> > & Ndir, int checkForcedYield = 0) const;
	
    void Read(TPZStream& buf, void* context) override {
        
    }
    void Write(TPZStream &buf, int withclassid) const override{
        
    }


	//////////////////CheckConv related methods/////////////////////
public:
	
	/** @brief Number of types of residuals */
	int NumCases()
	{
		return 2;
	}
	
	static TPZTensor<REAL> gRefTension;
	/** @brief LoadState will keep a given state as static variable of the class */
	void LoadState(TPZFMatrix<REAL> &state)
	{
		int i;
		for(i=0; i<6; i++) gRefTension.fData[i] = state(i,0);
	}
	
	void ComputeTangent(TPZFMatrix<REAL> &tangent, TPZVec<REAL> &, int icase)
	{
		TPZTensor<REAL> gradtheta;
		REAL theta,A(1.e7);
		TPZVec<REAL> phivec(1,0.);
		TPZVec<TPZTensor<REAL> > ndir(1);
		int i;
		switch(icase)
		{
			case 0:
				GradTheta(gRefTension,theta,gradtheta);
				tangent.Redim(1,6);
				for(i=0; i<6; i++) tangent(0,i) = gradtheta.fData[i];
				break;
			case 1:
				N(gRefTension,A,ndir,0);
				tangent.Redim(1,6);
				for(i=0; i<6; i++) tangent(0,i) = ndir[0].fData[i];
		}
	}
	
	void Residual(TPZFMatrix<REAL> &res,int icase)
	{
		TPZTensor<REAL> gradtheta;
		REAL theta,A(1.e7);
		TPZVec<REAL> phivec(1,0.);
		switch(icase)
		{
			case 0:
				GradTheta(gRefTension,theta,gradtheta);
				res.Redim(1,1);
				res(0,0) = theta;
				break;
			case 1:
				Compute(gRefTension,A,phivec,0);
				res.Redim(1,1);
				res(0,0) = phivec[0];
				break;
		}
		
	}

    virtual int GetNYield() const  override {
        return as_integer(NYield);
    }
};



template < class T>
void TPZYCTrescaRegularized::Compute(const TPZTensor<T> & sigma, const T & A,TPZVec<T> &result, int checkForcedYield) const
{
#ifdef PZ_LOG
	TPZLogger logger("plasticity.yctresca");
#endif
	
	//  result[0] = sqrt(sigma.J2()) - A;
	
	//  return;
	const REAL tol = 1.e-4;
	T invangle = InverseAngle(sigma);
	TPZTensor <T> s;
	sigma.S(s);
	REAL aux = (1. - tol);
	if (fabs(TPZExtractVal::val(invangle)) < aux)
	{
#ifdef PZ_LOG
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__ << " O angulo nao corresponde " << invangle;
		LOGPZ_DEBUG(logger,sout.str().c_str());
#endif
		TPZYCTresca::Compute(sigma,A,result,checkForcedYield);
	} else
	{
		if (fabs(TPZExtractVal::val(invangle)) > 1.) invangle = 1.;
		REAL asineps = asin(1. - tol);
		T tayext = ((T)fabs(TPZExtractVal::val(invangle)) - T(1.-tol))*T(1./(sqrt(1.- (1. - tol)*(1. - tol))));
		invangle = (tayext + T(asineps)) / 3.;
		//invangle = alpha * (1./3.) * (asin(1. - 1.e-6) + ( fabs(alpha)-(1.-1.e-6))*(1./(sqrt(1.-(1.-1.e-6)*(1.-1.e-6)))));
#ifdef PZ_LOG
		std::stringstream sout;
		T small = cos(invangle)*2.-T(sqrt(3.));
		
		sout << __PRETTY_FUNCTION__ << "2.*cos(invangle)-sqrt(3.) " << small;
		sout << " invangle " << invangle;
		LOGPZ_DEBUG(logger,sout.str().c_str());
#endif
		result[0] = sqrt(sigma.J2())*cos(invangle)*2. - A;
		//result[0] = sqrt(3.)* sqrt(sigma.J2()) - A;
	}
}

template <class T>
void TPZYCTrescaRegularized::GradTheta(const TPZTensor<T> & sigma,T & theta, TPZTensor<T> & gradtheta) const
{
#ifdef PZ_LOG
	TPZLogger logger("plasticity.yctresca");
#endif
	
	const REAL tol = 1.e-4;
	
	T invangle = InverseAngle(sigma);
	
	{
		/*    std::stringstream sout;
		 sout << "GradTheta:Invangle " << invangle;
		 LOGPZ_DEBUG(logger,sout.str().c_str());*/
	}
	
	if (fabs(TPZExtractVal::val(invangle)) < (1. - 1.e-6))
	{
#ifdef PZ_LOG
		{
			std::stringstream sout;
			sout << "calling GradTheta from father... invangle " << invangle;
			LOGPZ_DEBUG(logger,sout.str().c_str());
		}
#endif
		
		TPZYCTresca::GradTheta(sigma,theta,gradtheta);
#ifdef PZ_LOG
		{
			std::stringstream sout;
			sout << "sigma "  << sigma << std::endl;
			sout << "theta " << theta << std::endl;
			sout << "gratheta " << gradtheta << std::endl;
			LOGPZ_DEBUG(logger,sout.str().c_str());
		}
#endif
		
		
	} else
	{
		//double alpha;
		//if (fabs(invangle) > 1.)
		//     alpha = (invangle > 0.) ? 1. : -1.;
		//     theta = (alpha/3.)* (asin(1. - 1.e-6) + (fabs(alpha) - (1.-1.e-6))*(1./(sqrt(1.-(1.-1.e-6)*(1.-1.e-6)))));
		//     GradInverseAngle(sigma,gradtheta);
		//     T derivasin = 1./sqrt(1.-(1. - 1.e-6)*(1. - 1.e-6));
		//     gradtheta.Multiply(derivasin,1./3.);
		REAL alpha = (TPZExtractVal::val(invangle) > 0.) ? 1. : -1.;
		if (fabs(TPZExtractVal::val(invangle) ) > 1.) invangle = 1.;
		REAL asineps = asin(1. - tol);
		T tayext = ((T)fabs(TPZExtractVal::val(invangle)) - T(1.- tol))*T(1./(sqrt(1.- (1. - tol)*(1. - tol))));
		theta = (tayext + T(asineps)) * T(alpha / 3.);
		GradInverseAngle(sigma,gradtheta);
		REAL derivasin = 1./sqrt(1.-(1. - tol)*(1. - tol));
		gradtheta.Multiply(T(derivasin),T(1./3.));
		{
			/*      std::stringstream sout;
			 sout << "Invangle " << invangle << std::endl;
			 sout << "alpha    " << alpha << std::endl;
			 sout << "asineps  " << asineps << std::endl;
			 sout << "tayext   " << tayext << std::endl;
			 sout << "theta    " << theta << std::endl;
			 sout << "gradtheta" << gradtheta << std::endl;
			 sout << "derivsin " << derivasin << std::endl;
			 sout << "gradtheta" << gradtheta << std::endl;
			 LOGPZ_DEBUG(logger,sout.str().c_str());*/
		}
	}
}

template <class T> 
void TPZYCTrescaRegularized::N(const TPZTensor<T> & sigma,const T & A,  TPZVec<TPZTensor<T> > & Ndir, int checkForcedYield) const
{
#ifdef PZ_LOG
	TPZLogger logger("plasticity.yctresca");
#endif
	
	//   sigma.dJ2(Ndir[0]);
	//    Ndir[0].Multiply(sqrt(3.)*0.5*(1./sqrt(sigma.J2())),1.);
	//  return;
	
	
	
	/*  T invangle = InverseAngle(sigma);
	 T theta = (1./3.)*asin(invangle);
	 TPZTensor<T> gradinvangle;
	 GradInverseAngle(sigma,gradinvangle);
	 T derivasin = 1./sqrt(1.-invangle*invangle);*/
	//  T invangle;
	T theta;
	TPZTensor<T> gradtheta;
	//T derivasin;
	//GradASinInvAngle (sigma,invangle,theta,gradinvangle,derivasin);
	{
		/*    std::stringstream sout;
		 sout << "A  " << A << "    N: sigma " << sigma;
		 LOGPZ_DEBUG(logger,sout.str().c_str());*/
	}
	
	GradTheta<T> (sigma,theta,gradtheta);
	
	{
		/*    std::stringstream sout;
		 sout << "N: theta " << theta;
		 sout << "gradtheta " << gradtheta;
		 LOGPZ_DEBUG(logger,sout.str().c_str());*/
	}
	
	T sqj2 = sqrt(sigma.J2());
	TPZTensor<T> dj2;
	sigma.dJ2(dj2);
	
	TPZTensor<T> resultdtheta(gradtheta);
	resultdtheta.Multiply(sqj2*sin(theta),-2.);
	TPZTensor<T> resultdj2(dj2);
	resultdj2.Multiply(cos(theta)/sqj2,1.);
	TPZTensor<T> result(resultdj2);
	//   result.Add(resultdtheta,1.);
	
	//  resultdj2.Multiply(1./sqj2,sqrt(3.)/2.);
	
	
	{
		/*    std::stringstream sout;
		 sout << "sqj2 " << sqj2;
		 sout << "   theta" << theta;
		 T coco = sin (theta);
		 sout << "  sintheta " << coco;
		 sout << "dj2" << dj2;
		 sout << "N: resultdtheta " << theta;
		 LOGPZ_DEBUG(logger,sout.str().c_str());*/
	}
	{ 
		/*    std::stringstream sout;
		 sout << "gradtheta " << gradtheta;
		 sout << "Result N = " << result;
		 LOGPZ_DEBUG(logger,sout.str().c_str());*/
	}
	
	Ndir[0] = result;
	//  Ndir[0] = resultdj2;
}




#endif
