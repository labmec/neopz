//$Id: TPZYCDruckerPrager.h,v 1.5 2009-12-14 21:59:52 erick Exp $
#ifndef TPZYCDRUCKERPRAGER_H
#define TPZYCDRUCKERPRAGER_H

#include "TPZTensor.h"
#include "pzfmatrix.h"

#include "pzlog.h"
#include "pzsave.h"


#ifdef LOG4CXX // LOG4CXX may be defined alone or with LOG4CXX_PLASTICITY. The latter shall not be used alone.
#include <log4cxx/logger.h>
#include <log4cxx/basicconfigurator.h>
#include <log4cxx/propertyconfigurator.h>
#endif

#ifdef LOG4CXX_PLASTICITY
static LoggerPtr loggerDP(Logger::getLogger("plasticity.DruckerPrager"));
#endif




/**
 Implementa  a plastificacao do criterio de Von Mises
 */
class TPZYCDruckerPrager: public TPZSaveable {
    
	
public:
	
	
    TPZYCDruckerPrager():fKsi(0.),fEta(0.){}
	
    TPZYCDruckerPrager(const TPZYCDruckerPrager & source)
    {
		fKsi  = source.fKsi;
   		fEta  = source.fEta;
    }
	
	enum {NYield = 1};
	
    const char * Name() const
    {
		return "TPZYCDruckerPrager";	
    }

    void Print(std::ostream & out) const
    {
		out << Name();
    }
	
	/**
	 * Setup of material parameters
	 * @param [in] phi Mohr Coulomb's internal friction angle
	 * @param [in] innerMCFit If one, Drucker Prager model is inscribed in a referred Mohr Coulomb envelope. If zero, circumscribed.
	 * VERY IMPORTANT!! The ThermoForceA parameters should be set as:
	 *      fk: Herdening slope for the cohesion
	 *      fYield0 : equivalent Mohr Coulomb cohesion C
	 */
	void SetUp(const REAL & phi, const int innerMCFit) 
	{
		if(innerMCFit == 0)
		{
			fKsi = 6. * cos(phi)/(sqrt(3.) * (3.+sin(phi)));//INNER
			fEta = 6. * sin(phi)/(sqrt(3.)*(3.+sin(phi)));
			return;	
		}
			fKsi = 6. * cos(phi)/(sqrt(3.) * (3.-sin(phi)));//OUTER
			fEta = 6. * sin(phi)/(sqrt(3.)*(3.-sin(phi)));
	}
	
	int GetForceYield()
	{
		return 0; //nothing to be done in this yield criterium
	}
	
	void SetForceYield(const int forceYield)
	{
		//nothing to be done in this yield criterium
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
    void Compute(const TPZTensor<T> & sigma, const T & A, TPZVec<T> &res, int checkForcedYield = 0)const;
    
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
	
	virtual void Write(TPZStream &buf, int withclassid);
	virtual void Read(TPZStream &buf, void *context);
	virtual int ClassId() const;


public:
	REAL fPhi;

	
	
	
	REAL fKsi, fEta;

	//////////////////CheckConv related methods/////////////////////
};

template <class T>
void TPZYCDruckerPrager::Compute(const TPZTensor<T> & sigma,const T & A, TPZVec<T> &res, int checkForcedYield) const
{	    
		T I1,J2,p;
		I1 = sigma.I1();
	    J2 = sigma.J2();
		p = I1 / T(3.);
		
	    res[0] =  sqrt( J2 ) + p * T(fEta) - A * T(fKsi);
	
#ifdef LOG4CXX_PLASTICITY
	{
		std::stringstream sout;
		sout << "\nJ2 = " << J2;
		sout << "\nI1 = " << I1;	
		sout << "\nPHI = " << res[0];
		sout << " \nA = " << A;
		LOGPZ_INFO(loggerDP,sout.str());
	}
#endif

}

template <class T> 
void TPZYCDruckerPrager::N(const TPZTensor<T> & sigma,const T & A, TPZVec<TPZTensor<T> > & Ndir, int checkForcedYield) const
{
	//Deviatoric part
	T J2 = sigma.J2();
	TPZTensor<T> s;		
	sigma.S(s);
	s *= T(0.5) / sqrt(J2);
	//Hydrostatic part
	T EtaOver3 = T(fEta/3.);
	s.XX() += EtaOver3;
	s.YY() += EtaOver3;
	s.ZZ() += EtaOver3;
	Ndir[0] = s;
	
}

template <class T> 
void TPZYCDruckerPrager::H(const TPZTensor<T> & sigma,const T & A, TPZVec<T> & h, int checkForcedYield) const
{
	//H=-dphi/dA
    h[0] = fKsi;
}


inline void TPZYCDruckerPrager::Write(TPZStream &buf, int withclassid)
{
}
inline void TPZYCDruckerPrager::Read(TPZStream &buf, void *context)
{
}
inline int TPZYCDruckerPrager::ClassId() const
{
	return 888888;
}
inline template class TPZRestoreClass<TPZYCDruckerPrager, 888888> ;

#endif//TPZYDruckerPrager
