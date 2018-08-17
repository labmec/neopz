/**
 * @file
 */

#ifndef TPZYCDRUCKERPRAGER_H
#define TPZYCDRUCKERPRAGER_H

#include "TPZTensor.h"
#include "pzfmatrix.h"

#include "pzlog.h"
#include "TPZSavable.h"
#include "TPZPlasticCriterion.h"

#ifdef LOG4CXX // LOG4CXX may be defined alone or with LOG4CXX_PLASTICITY. The latter shall not be used alone.
#include <log4cxx/logger.h>
#include <log4cxx/basicconfigurator.h>
#include <log4cxx/propertyconfigurator.h>
#endif

#ifdef LOG4CXX_PLASTICITY
static LoggerPtr loggerDP(Logger::getLogger("plasticity.DruckerPrager"));
#endif

/**
 * @brief Implementa  a plastificacao do criterio de Von Mises
 */
class TPZYCDruckerPrager: public TPZPlasticCriterion {    
	
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

    void Print(std::ostream & out) const override
    {
		out << Name();
    }
	
	/**
	 * @brief Setup of material parameters
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
	 * @brief Checks if the proposed yield state leads to post-peak material behaviour. If so, the material
	 * is forced to behave in post-peak in order to avoid equation switching during Newton's method
	 * in the PlasticLoop routines.
	 * @param[in] sigma stress state
	 * @param[in] A Thermo Force
	 */
	void SetYieldStatusMode(const TPZTensor<REAL> & sigma, const REAL & A)
	{
		// nothing to be done in this yield criterium
	}
	
    /**
	 @brief Calculo do criterio de plastificacao 
	 @param[in] sigma tensao atual
	 @param[in] A forca thermodinamica atual
	 @param res
	 @param[in] checkForcedYield indicates wether to force post-peak failure behavior
	 */  
    template < class T>
    void Compute(const TPZTensor<T> & sigma, const T & A, TPZVec<T> &res, int checkForcedYield = 0)const;
    
    /**
	 @brief Derivada da funcao de plastificacao
	 @param[in] sigma tensao atual
	 @param[in] A forca termodinamica atual
	 @param[out] Ndir Derivada com respeito a tensao
	 @param[in] checkForcedYield indicates wether to force post-peak failure behavior
	 */
    template <class T> 
    void N(const TPZTensor<T> & sigma,const T & A,  TPZVec<TPZTensor<T> > & Ndir, int checkForcedYield = 0) const;
	
    /**
	 @brief Derivada da funcao de plastificacao com respeito a forca termodinamica
	 @param[in] sigma tensao atual
	 @param[in] A forca termodinamica atual
	 @param[out] h Derivada com respeito a forca termodinamica
	 @param[in] checkForcedYield indicates wether to force post-peak failure behavior
	 */
    template <class T> 
    void H(const TPZTensor<T> & sigma,const T & A,  TPZVec<T> & h, int checkForcedYield = 0) const;
	
    /**
     * Multiplicador para o caso onde utilizamos uma variavel de dano modificada
     */
    template <class T>
    void AlphaMultiplier(const T &A, T &multiplier) const
    {
        multiplier = T(1.);
    }
    
    
	virtual void Write(TPZStream &buf, int withclassid) const override;
	virtual void Read(TPZStream &buf, void *context) override;
	public:
virtual int ClassId() const override;

    void YieldFunction(const TPZVec<STATE>& sigma, STATE kprev, TPZVec<STATE>& yield) const override {
        TPZTensor<STATE> sigmaTensor;
        sigmaTensor.XX() = sigma[0];
        sigmaTensor.YY() = sigma[1];
        sigmaTensor.ZZ() = sigma[2];
        Compute(sigmaTensor, kprev, yield, 0);
    }

    virtual int GetNYield() const override {
        return as_integer(NYield);
    }


public:
	REAL fPhi;

	REAL fKsi, fEta;

	//////////////////CheckConv related methods/////////////////////
    
    
    /** @brief Number of types of residuals */
    int NumCases()
    {
        return 2;
    }
    TPZTensor<REAL> gRefTension;
    
    /** @brief LoadState will keep a given state as static variable of the class */
    void LoadState(TPZFMatrix<REAL> &state)
    {
        int i;
        for(i=0; i<6; i++) gRefTension.fData[i] = state(i,0);
    }
    void ComputeTangent(TPZFMatrix<REAL> &tangent, TPZVec<REAL> &coefs, int icase)
    {
        switch(icase)
        {
            case 0:
            {
                TPZVec<TPZTensor<REAL> > Ndir(1);
                REAL yield = 1.e6;
                this->N<REAL>(gRefTension,yield, Ndir, 0);
                //this->SetUp(20.,0);
                tangent.Redim(1,6);
                for(int i=0; i<6; i++)
                {
                    tangent(0,i) = Ndir[0].fData[i];
                }
                break;
            }
            case 1:
            {
                TPZTensor<REAL> dj2;
                gRefTension.dJ2(dj2);
                tangent.Redim(1,6);
                for(int i=0; i<6; i++)
                {
                    tangent(0,i) = dj2.fData[i];
                }
                break;
                
            }
                
        }
    }
    
    void Residual(TPZFMatrix<REAL> &res,int icase)
    {
        
        res.Redim(1,1);
        
        switch(icase)
        {
                
            case 0:
            {
                TPZVec<REAL> phi(1);
                REAL yield = 1.e6;
                //this->SetUp(20.,0);
                this->Compute(gRefTension,yield,phi,0);
                res(0,0) = phi[0];
                break;
            }
            case 1:
            {
             
                REAL j2 = gRefTension.J2();
                res(0,0)=j2;
            }
        }
        
    }
    
    

    
public:

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
	TPZTensor<T> s,dj2;


    
//SOUZA NETO CHECK CONV NAO DA 2 com S, mas com dJ2 BATE!
//    sigma.S(s);
//	s *= T(0.5) / sqrt(J2);
//Hydrostatic part
//	T EtaOver3 = T(fEta/3.);
//	s.XX() += EtaOver3;
//	s.YY() += EtaOver3;
//	s.ZZ() += EtaOver3;
//	Ndir[0] = s;
    
    sigma.dJ2(dj2);
    dj2*=T(0.5)/sqrt(J2);
    T EtaOver3 = T(fEta/3.);
    dj2.XX() += EtaOver3;
    dj2.YY() += EtaOver3;
    dj2.ZZ() += EtaOver3;
    Ndir[0]=dj2;
	
}

template <class T> 
void TPZYCDruckerPrager::H(const TPZTensor<T> & sigma,const T & A, TPZVec<T> & h, int checkForcedYield) const
{
	//H=-dphi/dA
    h[0] = fKsi;
}


inline void TPZYCDruckerPrager::Write(TPZStream &buf, int withclassid = 0) const
{
}
inline void TPZYCDruckerPrager::Read(TPZStream &buf, void *context = 0)
{
}

#endif//TPZYDruckerPrager
