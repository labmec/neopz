/**
 * @file
 */

#ifndef TPZYCLADEKIM_H
#define TPZYCLADEKIM_H

#include "TPZTensor.h"
#include "pzfmatrix.h"
#include "pzlog.h"

#ifndef CHECKCONV
#define CHECKCONV
#include "checkconv.h"
#include "fadType.h"
#include "TPZPlasticCriterion.h"
#endif

#ifdef PZ_LOG
static TPZLogger loggerYCLadeKim("plasticity.LadeKim");
#endif

/**
 * @brief Implementa as funções de potencial plástico e yield criterium do 
 * modelo constitutivo de Lade Kim para solo e rochas brandas
 */
class TPZYCLadeKim : public TPZPlasticCriterion {
    
public:

  enum {NYield = 1};
  
    virtual int ClassId() const override;

    TPZYCLadeKim():fKsi1(0.),fh(0.),m_hardening(0.),fKsi2(0.),fMu(0.),fNeta1(0.),fm(0.),fPa(0.),fForceYield(0){ }
	
    TPZYCLadeKim(const TPZYCLadeKim & source)
    {
		fKsi1  = source.fKsi1;
   		fh     = source.fh;
		m_hardening = source.m_hardening;
		fKsi2  = source.fKsi2;
		fMu    = source.fMu;
		fNeta1 = source.fNeta1;
		fm     = source.fm;
		fPa    = source.fPa;
		fForceYield = source.fForceYield;
    }

    TPZYCLadeKim & operator=(const TPZYCLadeKim & source)
    {
		fKsi1  = source.fKsi1;
   		fh     = source.fh;
		m_hardening = source.m_hardening;
		fKsi2  = source.fKsi2;
		fMu    = source.fMu;
		fNeta1 = source.fNeta1;
		fm     = source.fm;
		fPa    = source.fPa;
		fForceYield = source.fForceYield;
		return *this;
    }
    
    void Write(TPZStream& buf, int withclassid) const override {
        buf.Write(&fKsi1);
        buf.Write(&fh);
        buf.Write(&m_hardening);
        buf.Write(&fKsi2);
        buf.Write(&fMu);
        buf.Write(&fNeta1);
        buf.Write(&fm);
        buf.Write(&fPa);
        buf.Write(&fForceYield);
    }

    void Read(TPZStream& buf, void* context) override {
        buf.Read(&fKsi1);
        buf.Read(&fh);
        buf.Read(&m_hardening);
        buf.Read(&fKsi2);
        buf.Read(&fMu);
        buf.Read(&fNeta1);
        buf.Read(&fm);
        buf.Read(&fPa);
        buf.Read(&fForceYield);
    }
    
	const char * Name() const
    {
	   return "TPZYCLadeKim";	
    }
	
    void Print(std::ostream & out) const override
    {
		out << "\n" << this->Name();
		out << "\n fKsi1  = " << fKsi1;
		out << "\n fh     = " << fh;
		out << "\n m_hardening = " << m_hardening;
		out << "\n fKsi2  = " << fKsi2;
		out << "\n fMu    = " << fMu;
		out << "\n fNeta1 = " << fNeta1;
		out << "\n fm     = " << fm;
		out << "\n fPa    = " << fPa;
		out << "\n fForceYield = " << fForceYield;
    }
	
	int GetForceYield()
	{
		return fForceYield;
	}
	
	void SetForceYield(const int forceYield)
	{
		fForceYield = forceYield;
	}
	
    /**
	 * @brief Calculo do criterio de plastificacao 
	 * @param[in] sigma tensao atual
	 * @param res vector
	 * @param[in] A forca thermodinamica atual
	 * @param[in] checkForcedYield indicates wether to force post-peak failure behavior
	 */  
    template < class T>
    void Compute(const TPZTensor<T> & sigma, const T & A, TPZVec<T> &res, int checkForcedYield) const;

    /**
	 * @brief Calculo da função de potencial plástico
	 * @param[in] sigma tensao atual
	 * @param[in] A forca thermodinamica atual
	 * @param[out] PlasticPot Derivada com respeito a tensao
	 * @param[in] checkForcedYield indicates wether to force post-peak failure behavior
	 */  
    template < class T>
    void ComputePlasticPotential(const TPZTensor<T> & sigma, const T & A, T & PlasticPot, int checkForcedYield) const;

    /**
	 * @brief Derivada da derivada da funcao de potencial plastico (direção de plastificação)
	 * @param[in] sigma tensao atual
	 * @param[in] A forca termodinamica atual
	 * @param[out] Ndir Derivada com respeito a tensao
	 * @param[in] checkForcedYield indicates wether to force post-peak failure behavior
	 */
    template <class T> 
    void N(const TPZTensor<T> & sigma,const T & A,  TPZVec<TPZTensor<T> > & Ndir, int checkForcedYield) const;

    /**
	 * @brief Derivada da funcao de plastificacao com respeito a forca termodinamica
	 * @param[in] sigma tensao atual
	 * @param[in] A forca termodinamica atual
	 * @param[out] h Derivada com respeito a forca termodinamica
	 * @param[in] checkForcedYield indicates wether to force post-peak failure behavior
	 */
    template <class T> 
    void H(const TPZTensor<T> & sigma,const T & A,  TPZVec<T> & h, int checkForcedYield) const;

    inline void SetUp(REAL Ksi1, REAL Ksi2, REAL h, REAL Alpha, REAL Mu, REAL neta1, REAL m, REAL pa)
    {
       fKsi1  = Ksi1;
       fKsi2  = Ksi2;
       fh     = h;
       m_hardening = Alpha;
       fMu    = Mu;
       fNeta1 = neta1;
       fm     = m;
       fPa    = pa;
	
    }

    /**
     * Multiplicador para o caso onde utilizamos uma variavel de dano modificada
     */
    template <class T>
    void AlphaMultiplier(const T &A, T &multiplier) const
    {
        multiplier = T(1.);
    }

	/**
	 * @brief Checks if the proposed yield state leads to post-peak material behaviour. If so, the material
	 * is forced to behave in post-peak in order to avoid equation switching during Newton's method
	 * in the PlasticLoop routines.
	 * @param[in] sigma stress state
	 * @param[in] A Thermo Force
	 */
	void SetYieldStatusMode(const TPZTensor<REAL> & sigma, const REAL & A);
        
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

   /**
	* @brief Parameter related to the YC and Plastic Potential
	* Ksi1 models the shape of the yield funcition at the deviatoric planes by defining the relative influence of invariants I3 and I2.
	* The I3 term leads to a triangular shaped YC at the deviatoric planes while the I2 term leads to a circular shaped YC.
	* OBS: Ksi1 may be related to the m parameter (curvature of meridians for failure surface)
	* through the expression: ksi1 = 0.00155 * m ^ -1.27
	*/
	REAL fKsi1;

   /**
	* @brief Parameter related to the YC \n
	* The h constant models the curvature of the YC meridians, i.e., how the meridians
	* vary along with the level of hydrostatic stress (I1/3).
	*/
	REAL fh;

   /**
     @brief Parameter related to the YC \n
     Alpha models how the exponent q varies according to the proximity of the stress state
     to the failure surface (surface at which the material starts to soften).
     In the bibliography, q is defined as a function of the stress level only.
     During isotropic compression, q evaluates to null and the material is allowed
     to plastify and harden unlimitedly;
     Stress paths that lead to deviatoric stresses correspond to values of q between 0 and 1,
     assuming 1 when the material reaches the failure surface.
     Phisically one may say that the parameter Alpha (and therefore q too) is responsible
     for the material behavior of increasing rates of hardening with increasing deviatoric
     stresses.
     Notice that this alpha is not related to the Plastic Damage variable.
   */
   REAL m_hardening;

   /**
     @brief Parameter related to the Plastic Potential \n
     Ksi2 controls the intersection of the Plastic Potential with the Hydrostatic axis.
   */
   REAL fKsi2;

   /**
     @brief Parameter related to the Plastic Potential \n
     Mu defines the curvature of the meridians.
   */
   REAL fMu;

   /**
     @brief Parameter related to the Failure Surface \n
     Neta1 is the value of the failure Surface for the material when it fails.
   */
   REAL fNeta1;

   /**
     @brief Parameter related to the Failure Surface \n
     m models the curvature of the meridians of the Failure Surface.
   */
   REAL fm;

   /** @brief Atmospheric pressure to input/remove dimensional effects */
   REAL fPa;
	
   /** @brief Post Peak material behavior */
   int fForceYield;

public:
//////////////////CheckConv related methods/////////////////////

    /** @brief number of types of residuals */
    inline int NumCases();

    static TPZTensor<REAL> gRefTension;
    /** @brief LoadState will keep a given state as static variable of the class */
    inline void LoadState(TPZFMatrix<REAL> &state);

    inline void ComputeTangent(TPZFMatrix<REAL> &tangent, TPZVec<REAL> &, int icase);

    inline void Residual(TPZFMatrix<REAL> &res,int icase);

    static void CheckConv()
{
   const int nVars = 6;

   /** @brief Creating the LadeNelsonElasticResponse obejct */
   TPZYCLadeKim YCLadeKim;
   // setup with data from Sacramento River Sand
   // Lade, Poul V., Kim, Moon K. Single Hardening Constitutive Model
   // for soil, Rock and Concrete. Int. J. Solids Structures, vol32 
   // No14 pp 1963-1978,1995
   REAL m     = 0.093,
        Ksi2  = -3.09,
        h     = 0.765,
        Alpha = 0.229,
        Mu    = 2.00,
        neta1 = 80.,
        pa    = 14.7;
   REAL Ksi1  = 0.00155 * pow(m, -1.27);

   YCLadeKim.SetUp(Ksi1, Ksi2, h, Alpha, Mu, neta1, m, pa);

   TPZFNMatrix<nVars> sigma(nVars,1), Range(nVars,1);
   sigma(_XX_,0) = 0.17*100.;
   sigma(_YY_,0) = 0.13*100.;
   sigma(_ZZ_,0) = 0.11*100.;
   sigma(_XY_,0) = 0.7 *100.;
   sigma(_XZ_,0) = 0.5 *100.;
   sigma(_YZ_,0) = 0.3 *100.;

   Range = sigma * (1./19.);
   TPZVec< REAL > Coefs(1,1.);
   CheckConvergence(YCLadeKim, sigma, Range, Coefs);
}
//////////////////CheckConv related methods/////////////////////
    
public:
 

};


template <class T>
inline void TPZYCLadeKim::Compute(const TPZTensor<T> & sigma,const T & A, TPZVec<T> &res, int checkForcedYield) const
{
	bool negI1 = false, output;
    T S, q;
    T I1 = sigma.I1();
	if( fabs(TPZExtractVal::val(I1) ) < 1.e-24) I1 /*+*/= 1.e-24; // avoiding Log(0)
	if( TPZExtractVal::val(I1) < 0.)
	{
		I1 = -I1; // avoiding Log(Neg)
		negI1 = true;
	}
    T I2 = sigma.I2();
	if( fabs(TPZExtractVal::val(I2) ) < 1.e-24) I2 /*+*/= 1.e-24; // avoiding division by zero
    T I3 = sigma.I3();
	if( fabs(TPZExtractVal::val(I3) ) < 1.e-24) I3 /*+*/= 1.e-24; // avoiding division by zero
	//if( TPZExtractVal::val(I3) < 0 ) I3 -= 2* TPZExtractVal::val(I3);

    T I13_I3 = I1 * I1 * I1 / I3;
    T I12_I2 = I1 * I1 / I2;

    //S = T( 1./fNeta1)  * ( I13_I3 - T(27.) ) * pow( I1/T(fPa), T(fm) );
	if(!negI1)
	{
		S = T( 1./fNeta1)  * ( I13_I3 - T(27.) ) * exp( T(fm) * log( I1/T(fPa) ) );
	}else
	{
		S = T(1.);
	}
	REAL S_real = TPZExtractVal::val(S);
	
	if(checkForcedYield && fForceYield)
	{
		#ifdef PZ_LOG
    	{
			std::stringstream sout;
			sout << "*** Compute *** Imposing S = 1.0 (when it's " << S_real 
				 << ") because fForcedYield was set to TRUE within this PlasticLoop";
			LOGPZ_INFO(loggerYCLadeKim,sout.str().c_str());
	    }
		#endif
		S = T(1.);
		q = T(1.);
	}else
	{//checkForcedYield && fForceYield == FALSE
		output = false;
	
		if(S_real > 1.)
		{
			output = true;
			S = T(1.);
			// When S > 1 it means that the material has already failed. This implementation
			// does not consider material softening like the original LadeKim Proposal because
			// such material behaviour is known to lead to mesh-dependent solutions, which is
			// not really desired in a FEM method. In this implementation the forced S=1 causes
			// the material to behave as a perfect plastic material when at failure stress state.
		}
	
		if(S_real < 0.)
		{
			output = true;
			if(S_real > -0.001)
			{
				S = T(0.);
			}else
			{
				S = T(1.);
				// By analysing the expression for S it can be seen that S may be neagtive only in 
				// the case where I3<0. That implies the whole stress state is negative - tension.
				// The main Lade-Kim hypothesis is that the stress state is always positive, always
				// in compression. In such an inconsistence, S is forced to assume the value of 1
				// causing the material to behave as a perfect plastic material and forcing the stress
				// state moving back at the failure state.
			}
		}
	
		if(output)
		{
			#ifdef PZ_LOG
   		 	{
				std::stringstream sout;
				sout << "** Compute *** Forcing S = " << TPZExtractVal::val(S) << " when S = " << S_real
					 << ".\nI1 = " << I1 << "\nI2 = " << I2 << "\nI3 = " << I3 
					 << "\nsigma = " << sigma;
				LOGPZ_WARN(loggerYCLadeKim,sout.str().c_str());
		    }
			#endif
		}
			
		q = T(m_hardening) * S / (T(1.) - T(1.-m_hardening) * S);
	}
	
	REAL q_real = TPZExtractVal::val(q);
	output = false;
	
    if( q_real < -1.e-10 )
	{
		output = true;
		q = T(0);
	}
	if( q_real > 1)
	{
		output = true;
		q = T(1);
	}
		
    if( output )
	{
		#ifdef PZ_LOG
    	{
			std::stringstream sout;
			sout << "** Compute *** Forcing q = " << q << " when q = " << q_real
				 << ".\nI1 = " << I1 << "\nI2 = " << I2 << "\nI3 = " << I3 
				 << "\nsigma = " << sigma;
			LOGPZ_WARN(loggerYCLadeKim,sout.str().c_str());
	    }
		#endif
		q = T(1.);
	}
	
    //res[0] = ( T(fKsi1) * I13_I3 - I12_I2 ) * pow( I1/ T(fPa), T(fh) ) * T(exp( q )) - A;
	res[0] = ( T(fKsi1) * I13_I3 - I12_I2 ) * exp( T(fh) * log( I1/ T(fPa) ) ) * T(exp( q )) - A;
}


template < class T>
inline void TPZYCLadeKim::ComputePlasticPotential(const TPZTensor<T> & sigma, const T & A,  T & PlasticPot, int checkForcedYield) const
{
    T I1 = sigma.I1();
	if( fabs(TPZExtractVal::val(I1) ) < 1.e-24) I1 /*+*/= 1.e-24; // avoiding Log(0)
	if( TPZExtractVal::val(I1) < 0.) I1 = -I1; // avoiding Log(Neg)
    T I2 = sigma.I2();
	if( fabs(TPZExtractVal::val(I2) ) < 1.e-24) I2 /*+*/= 1.e-24; // avoiding division by zero
    T I3 = sigma.I3();
	if( fabs(TPZExtractVal::val(I3) ) < 1.e-24) I3 /*+*/= 1.e-24; // avoiding division by zero
	//if( TPZExtractVal::val(I3) < 0 ) I3 -= 2* TPZExtractVal::val(I3);
	
    T I1_2 = I1 * I1;
    // PlasticPot = ( I1 * I1_2 / I3 * T(fKsi1) - I1_2 / I2 + T(fKsi2) ) * pow( I1/T(fPa), fMu );
	PlasticPot = ( I1 * I1_2 / I3 * T(fKsi1) - I1_2 / I2 + T(fKsi2) ) *
		exp( T(fMu) * log(I1/T(fPa) ) );
}

template <class T> 
inline void TPZYCLadeKim::N(const TPZTensor<T> & sigma, const T & A, TPZVec<TPZTensor<T> > & Ndir, int checkForcedYield) const
{
    T I1 = sigma.I1();
	if( fabs(TPZExtractVal::val(I1) ) < 1.e-24) I1 /*+*/= 1.e-24; // avoiding division by zero
	if( TPZExtractVal::val(I1) < 0.) I1 = -I1; // avoiding Log(Neg)
    T I2 = sigma.I2();
	if( fabs(TPZExtractVal::val(I2) ) < 1.e-24) I2 /*+*/= 1.e-24; // avoiding division by zero
    T I3 = sigma.I3();
	if( fabs(TPZExtractVal::val(I3) ) < 1.e-24) I3 /*+*/= 1.e-24; // avoiding division by zero
	//if( TPZExtractVal::val(I3) < 0 ) I3 -= 2* TPZExtractVal::val(I3);

    // I1_I2   = I1/I2
    // I12_I22 = I1²/I2²
    T I1_I2 = I1 / I2;
    T I12_I22 = I1_I2 * I1_I2;

    // I12_I3  = I1²/I3
    // I13_I32 = I1³/I3²
    T I12_I3 = I1 * I1 / I3;
    T ksi1_I13_I32 = I12_I3 * I1 / I3 * T(fKsi1);

    //T I1_Mu = pow( I1 / T(fPa) , fMu );
	T I1_Mu = exp( T(fMu) * log(I1 / T(fPa) ));
    T Two_I1_Mu = I1_Mu * T(2.);

    T G = T( fKsi1 * (fMu + 3.) ) * I12_I3 - 
          T( fMu + 2. ) * I1_I2 + 
          T( fKsi2 * fMu ) / I1;

    Ndir[0].XX() = I1_Mu * ( G - ( sigma.YY() + sigma.ZZ() ) * I12_I22
                               - ( sigma.YY() * sigma.ZZ() - sigma.YZ()*sigma.YZ() ) * ksi1_I13_I32 );

    Ndir[0].YY() = I1_Mu * ( G - ( sigma.ZZ() + sigma.XX() ) * I12_I22
                               - ( sigma.ZZ() * sigma.XX() - sigma.XZ()*sigma.XZ() ) * ksi1_I13_I32 );

    Ndir[0].ZZ() = I1_Mu * ( G - ( sigma.XX() + sigma.YY() ) * I12_I22
                               - ( sigma.XX() * sigma.YY() - sigma.XY()*sigma.XY() ) * ksi1_I13_I32 );

    Ndir[0].YZ() = Two_I1_Mu * ( I12_I22 * sigma.YZ()
                                 - ( sigma.XY() * sigma.XZ() - sigma.XX()*sigma.YZ() ) * ksi1_I13_I32 );

    Ndir[0].XZ() = Two_I1_Mu * ( I12_I22 * sigma.XZ()
                                 - ( sigma.XY() * sigma.YZ() - sigma.YY()*sigma.XZ() ) * ksi1_I13_I32 );

    Ndir[0].XY() = Two_I1_Mu * ( I12_I22 * sigma.XY()
                                 - ( sigma.YZ() * sigma.XZ() - sigma.ZZ()*sigma.XY() ) * ksi1_I13_I32 );



}

template <class T> 
inline void TPZYCLadeKim::H(const TPZTensor<T> & sigma,const T & A, TPZVec<T> & h, int checkForcedYield) const
{
    T PlasticPot;

    ComputePlasticPotential(sigma, A, PlasticPot, checkForcedYield);

    h[0] = PlasticPot * T(fMu);
}

inline void TPZYCLadeKim::SetYieldStatusMode(const TPZTensor<REAL> & sigma, const REAL & A)
{
	fForceYield = 0; // unsetting forceYield

    REAL S;
    REAL I1 = sigma.I1();
	//if( fabs(I1) < 1.e-24) I1 /*+*/= 1.e-24; // avoiding Log(0)
	if( I1 < 0.)
	{
		fForceYield = 1;
		#ifdef PZ_LOG
    	{
			std::stringstream sout;
			sout << "<<< SetYieldStatusMode *** Imposing fForceYield = TRUE because proposed S = " << S;
			LOGPZ_INFO(loggerYCLadeKim,sout.str().c_str());
	    }
		#endif
		return;
	}
    REAL I3 = sigma.I3();
	if( fabs(TPZExtractVal::val(I3) ) < 1.e-24) I3 /*+*/= 1.e-24; // avoiding division by zero

    REAL I13_I3 = I1 * I1 * I1 / I3;

    S = 1./fNeta1  * ( I13_I3 - 27. ) * pow( I1/fPa, fm );
	
	if(S < -1.e-10 || S > 1.)
	{
		fForceYield = 1;
		#ifdef PZ_LOG
    	{
			std::stringstream sout;
			sout << "<<< SetYieldStatusMode *** Imposing fForceYield = TRUE because proposed S = " << S;
			LOGPZ_INFO(loggerYCLadeKim,sout.str().c_str());
	    }
		#endif
		return;
	}
	
	#ifdef PZ_LOG
    {
		std::stringstream sout;
		sout << "<<< SetYieldStatusMode *** Leaving fForceYield = FALSE because proposed S = " << S;
		LOGPZ_INFO(loggerYCLadeKim,sout.str().c_str());
	}
	#endif

}

//////////////////CheckConv related methods/////////////////////

inline int TPZYCLadeKim::NumCases()
{
    return 2;
}

inline void TPZYCLadeKim::LoadState(TPZFMatrix<REAL> &state)
{
#ifdef PZ_LOG
    TPZLogger logger("plasticity.ycladekim");
#endif
  int i;
  for(i=0; i<6; i++) gRefTension[i] = state(i,0);
#ifdef PZ_LOG
  std::stringstream sout;
  sout << "Tension " << state;
  LOGPZ_DEBUG(logger,sout.str().c_str());
#endif
}

inline void TPZYCLadeKim::ComputeTangent(TPZFMatrix<REAL> &tangent, TPZVec<REAL> &, int icase)
{
#ifdef PZ_LOG
    TPZLogger logger("plasticity.ycladekim");
#endif

  const int nVars = 6;
  typedef TFad<nVars,REAL> TFAD;

  int i, j;
  TPZVec< TPZTensor < REAL > > N_Dir(1);
  REAL A; // not used, just created to satisfy interface
  TPZTensor < TFAD > Sigma_FAD;
  TFAD A_FAD;
  TPZVec< TPZTensor < TFAD > > N_Dir_FAD(1);

  switch(icase)
  {
    case 0:
      //Compute N
      tangent.Redim(1,nVars);
      N(gRefTension, A, N_Dir, 0);
      for(i=0; i<nVars; i++)
          tangent(0,i) = N_Dir[0][i];
      break;
    case 1:
      //Compute derivatives of N
      tangent.Redim(nVars,nVars);
      gRefTension.CopyTo(Sigma_FAD);
      for(i = 0;i < nVars; i++)
         Sigma_FAD[i].diff(i,nVars);
      N(Sigma_FAD, A_FAD, N_Dir_FAD, 0);
      for(i = 0; i < nVars; i++)
        for(j = 0; j < nVars; j++)
          tangent(i,j) = N_Dir_FAD[0][i].dx(j);
    break;

  }
#ifdef PZ_LOG
  std::stringstream sout;
  sout << "Matriz tangent " << tangent;
  LOGPZ_DEBUG(logger,sout.str().c_str());
#endif
}

inline void TPZYCLadeKim::Residual(TPZFMatrix<REAL> &res,int icase)
{
#ifdef PZ_LOG
    TPZLogger logger("plasticity.ycladekim");
#endif
  int i;
  const int nVars = 6;
//  typedef TFad<nVars,REAL> TFAD;

  REAL PlasticPot;
  REAL A; // not used, just created to satisfy interface
  TPZVec< TPZTensor<REAL> > N_Dir(1);
  switch(icase)
  {
    case 0:
      //Compute PlasticPotential
      res.Redim(1,1);
      ComputePlasticPotential(gRefTension, A, PlasticPot, 0);
      res(0,0) = PlasticPot;
    break;
    case 1:
      //Compute Ndir
      res.Redim(nVars,1);
      N(gRefTension, A, N_Dir, 0);
      for(i = 0; i < nVars; i++)
         res(i,0) = N_Dir[0][i];
    break;
  }
#ifdef PZ_LOG
  std::stringstream sout;
  sout << "Residual vector " << res;
  LOGPZ_DEBUG(logger,sout.str().c_str());
#endif
}


//////////////////CheckConv related methods/////////////////////


#endif //TPZYCLADEKIM_H
