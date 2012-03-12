// $Id: TPZYCSandlerDimaggio.h,v 1.11 2009-06-29 22:54:01 erick Exp $

#ifndef TPZYCSANDLERDIMAGGIO_H
#define TPZYCSANDLERDIMAGGIO_H

#include "TPZTensor.h"
#include "pzfmatrix.h"
#include "pzlog.h"

#ifndef CHECKCONV
#define CHECKCONV
#include "checkconv.h"
#endif

#include "fadType.h"

/**
Implementa as funções de potencial plástico e yield criterium do 
modelo constitutivo associativo de Sandler e Dimaggio (1971), desenvolvido
inicialmente para arenitos (Ranch McCormic Sand)
*/
class TPZYCSandlerDimaggio {
    

public:

  enum {NYield = 2};
	
    TPZYCSandlerDimaggio():fA(0.),fB(0.),fC(0.),fD(0.),fW(0.),fR(0.){ }
	
    TPZYCSandlerDimaggio(const TPZYCSandlerDimaggio & source)
    {
		fA = source.fA;
		fB = source.fB;
		fC = source.fC;
		fD = source.fD;
		fW = source.fW;
		fR = source.fR;
    }

    TPZYCSandlerDimaggio & operator=(const TPZYCSandlerDimaggio & source)
    {
		fA = source.fA;
		fB = source.fB;
		fC = source.fC;
		fD = source.fD;
		fW = source.fW;
		fR = source.fR;
		return *this;
    }

	const char * Name() const
    {
	   return "TPZYCSandlerDimaggio";	
    }
	
    void Print(std::ostream & out) const
    {
		out << "\n" << this->Name();
		out << "\n fA = " << fA;
		out << "\n fB = " << fB;
		out << "\n fC = " << fC;
		out << "\n fD = " << fD;
		out << "\n fR = " << fR;
		out << "\n fW = " << fW;
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
    void Compute(const TPZTensor<T> & sigma, const T & A, TPZVec<T> &res, int checkForcedYield) const;

    /**
    Derivada da derivada da funcao de potencial plastico (direção de plastificação)
    @param [in] sigma tensao atual
    @param [in] A forca termodinamica atual
    @param [out] Derivida com respeito a tensao
	@param [in] checkForcedYield indicates wether to force post-peak failure behavior
    */
    template <class T> 
    void N(const TPZTensor<T> & sigma,const T & A,  TPZVec<TPZTensor<T> > & Ndir, int checkForcedYield) const;

    /**
    Derivada da funcao de plastificacao com respeito a forca termodinamica
    @param [in] sigma tensao atual
    @param [in] A forca termodinamica atual
    @param [out] Derivida com respeito a forca termodinamica
	@param [in] checkForcedYield indicates wether to force post-peak failure behavior
    */
    template <class T> 
    void H(const TPZTensor<T> & sigma,const T & A,  TPZVec<T> & h, int checkForcedYield) const;

    inline void SetUp(REAL A, REAL B, REAL C, REAL D, REAL R, REAL W)
    {
	   fA = A;
	   fB = B;
	   fC = C;
	   fD = D;
	   fR = R;
	   fW = W;
    }

private:
   /**
    Solves for the invariant I1 value at the intersection
	of shear and hardening cap yield criteria. 
	The current value of L when entering the funcition is assumed
	to be the initial guess. 
	In this implementation L should be negative in compression.
   */
   template <class T>
   void SolveL(const T & X, T & L, REAL relTol = 1.e-6) const;

   /**
    Might be a reasonable initial guess for L when no better data is available. 
	In this implementation L should be negative in compression.
   */
   template <class T>
   void LInitialGuess(const T & X, T & L) const;
	
   /**
    Evaluates the F(L) grouping
   */
   template <class T>
   void ComputeF(const T & L, T & F) const;
	
   /**
    Evaluates the F(L) grouping total derivative
   */
   template <class T>
   void ComputedF(const T & L, T & dF) const;

   /**
    Evaluates X(EpsilonPvol), the value of the first invariant
	of an hidrostatic stress tensor at the cap yield surface.
	In this implementation X should negative in compression.
   */
   template <class T>
   void ComputeX(const T & A, T & X) const;
	
public:

		
   /**
     Parameter related to the YC and Plastic Potential
     A, B and C are constants in the modified Drucker-Prager
	 shear yield criteria.
   */
   REAL fA, fB, fC;

   /**
     Parameters related to the YC and Plastic Potential
	 The D and W parameters correlates the total plastic strain to the
	 hydrostatic loading level. It is thus related to the cap
	 hardening/softening.
   */
   REAL fD, fW;

   /**
     Parameter related to the YC and Plastic Potential
     half-axis ratio for the ellipsoidal hardening/softening
	 cap.
   */
   REAL fR;


public:
//////////////////CheckConv related methods/////////////////////

    /**
    number of types of residuals
    */
    inline int NumCases();

    static TPZTensor<REAL> gRefTension;
    /**
    LoadState will keep a given state as static variable of the class
    */
    inline void LoadState(TPZFMatrix &state);

    inline void ComputeTangent(TPZFMatrix &tangent, TPZVec<REAL> &, int icase);

    inline void Residual(TPZFMatrix &res,int icase);

    static void CheckConv();
	
//////////////////CheckConv related methods/////////////////////

//////////////////Internal routines verification/////////////////

	static void TestSolveL();
//////////////////Internal routines verification/////////////////

	static void McCormicRanchSand(TPZYCSandlerDimaggio & material);
};

template <class T>
inline void TPZYCSandlerDimaggio::Compute(const TPZTensor<T> & sigma,const T & A, TPZVec<T> &res, int checkForcedYield) const
{
	// the termoforce A in this case is assumed to be the
	// plastic volumetric strain itself. In fact it is not,
	// but the resultant derivatives are correct for practical purposes.
	
	// The following line evaluates L, the value of the first 
	// invariant of stresses at the intersection of the
	// shear and hardening cap yield criteria / plastic potential.
	// It is first evaluated as REAL type to avoid unnecessary
	// derivatives evaluation.
	
    #ifdef LOG4CXX_PLASTICITY
        {
          LoggerPtr logger(Logger::getLogger("plasticity.SandlerDimaggio"));
          std::stringstream sout;
          sout << ">>> TPZYCSandlerDimaggio::Compute *** - Plastic Potential / Yield - associative";
          LOGPZ_INFO(logger,sout.str().c_str());
        }
    #endif

	REAL L_REAL, X_REAL;
	ComputeX((REAL)shapeFAD::val(A), X_REAL);
	LInitialGuess(X_REAL, L_REAL);
	SolveL(X_REAL, L_REAL);
    
    T I1 = sigma.I1();
    T J2 = sigma.J2();
	
	// f1 - Modified Drucker-Prager as shear Yield Criterium
	T FI1;
	ComputeF(I1, FI1);
	if(fabs((REAL)shapeFAD::val(J2)) < 1.e-6)
	{
		res[0] = - FI1;	// avoiding nan derivatives
	}else{
		res[0] = sqrt(J2) - FI1;
	}
	
	// f2 - ellipsoidal hardening/softening cap
	T FL, L(L_REAL), X;
   	ComputeX(A, X);
	SolveL(X, L); // evaluating the derivatives of L
	ComputeF(L, FL);
		
	if(fabs( (double)shapeFAD::val(FL) ) < 0.00001)
	{
    #ifdef LOG4CXX_PLASTICITY
        {
          LoggerPtr logger(Logger::getLogger("plasticity.SandlerDimaggio"));
          std::stringstream sout;
          sout << "*** TPZYCSandlerDimaggio::ComputePlasticPotential ***";
          sout << "\nDivision by F=" << shapeFAD::val(L) << " at f2 - ellipsoidal hardening/softening cap";
          LOGPZ_WARN(logger,sout.str().c_str());
        }
    #endif
	}
	
	T Temp1( (L - I1)/(FL * T(fR) ) );
	Temp1 *= Temp1;
	T Temp2 = J2 / FL / FL;
		
	res[1] = Temp1 + Temp2 - T(1.);
	
	return;
}

template <class T> 
inline void TPZYCSandlerDimaggio::N(const TPZTensor<T> & sigma, const T & A, TPZVec<TPZTensor<T> > & Ndir, int checkForcedYield) const
{

	// the termoforce A in this case is assumed to be the
	// plastic volumetric strain itself. In fact it is not,
	// but the resultant derivatives are correct for practical purposes.
	
	//The following line evaluates L, the value of the first 
	// invariant of stresses at the intersection of the
	// shear and hardening cap yield criteria / plastic potential.
	// It is first evaluated as REAL type to avoid unnecessary
	// derivatives evaluation.
	
    #ifdef LOG4CXX_PLASTICITY
        {
          LoggerPtr logger(Logger::getLogger("plasticity.SandlerDimaggio"));
          std::stringstream sout;
          sout << ">>> TPZYCSandlerDimaggio::N *** - Plastification direction - associative";
          LOGPZ_INFO(logger,sout.str().c_str());
        }
    #endif
	
	REAL ResTol = 1.e-6;
	
	REAL L_REAL, X_REAL;
	ComputeX((double)shapeFAD::val(A), X_REAL);
	LInitialGuess(X_REAL, L_REAL);
	SolveL(X_REAL, L_REAL, ResTol);
    
    T I1 = sigma.I1();
    T J2 = sigma.J2();
    T SQRTJ2 = sqrt(J2);
	
	{
		//f1 - Modified Drucker-Prager as shear Yield Criterium / Plastic Potential

		T Temp1 = I1 * T(fB);
		Temp1 = exp( Temp1 ) * T (fB * fC);
		
		if((double)shapeFAD::val(SQRTJ2) < 1.e-6) // just for robustness. f1 shouldn't be reached when J2 = 0.
		{
			#ifdef LOG4CXX_PLASTICITY
            {
               LoggerPtr logger(Logger::getLogger("plasticity.SandlerDimaggio"));
               std::stringstream sout;
               sout << "*** TPZYCSandlerDimaggio::N *** - SQRT(J2) = " << shapeFAD::val(SQRTJ2) <<  " < 1.e-6 causes error in 0-th yield function. Imposing J2 = 1.e-6 instead";
               LOGPZ_WARN(logger,sout.str().c_str());
            }
            #endif
			SQRTJ2 = T(1.e-6);
		}
		
		Temp1 = Temp1 - I1 / SQRTJ2 / T(6.);
		T Temp2 = T(1.) / SQRTJ2;
		T Temp3 = Temp2 / T(2.);
		
		Ndir[0].XX() = Temp1 + sigma.XX() * Temp3;
		Ndir[0].YY() = Temp1 + sigma.YY() * Temp3;
		Ndir[0].ZZ() = Temp1 + sigma.ZZ() * Temp3;
	    Ndir[0].YZ() = sigma.YZ() * Temp2;
	    Ndir[0].XZ() = sigma.XZ() * Temp2;
	    Ndir[0].XY() = sigma.XY() * Temp2;
	}
	
	{//f2 - ellipsoidal hardening/softening cap

		T FL, X, L(L_REAL * 1.- ResTol); // guaranteeing that the function will be evaluated
	   	ComputeX(A, X);
		SolveL(X, L, ResTol); // evaluating the derivatives of L
		
		ComputeF(L, FL);
		T FL2 = FL * FL;
		T FL3 = FL2 / T(2.);
	
		T Temp = (I1-L)/ T(fR * fR) - I1 / T(6.);
		Temp = Temp / FL2 * T(2.);

			#ifdef LOG4CXX_PLASTICITY
            {
               LoggerPtr logger(Logger::getLogger("plasticity.SandlerDimaggio"));
               std::stringstream sout;
               sout << "*** TPZYCSandlerDimaggio::N *** X = " << X
					<< "\n L = " << L << " L_REAL = " << L_REAL
					<< "\n FL = " << FL
					<< "\n Temp = " << Temp;
               LOGPZ_DEBUG(logger,sout.str().c_str());
            }
            #endif
		
		Ndir[1].XX() = Temp + sigma.XX() / FL2;
		Ndir[1].YY() = Temp + sigma.YY() / FL2;
		Ndir[1].ZZ() = Temp + sigma.ZZ() / FL2;
        Ndir[1].YZ() = sigma.YZ() / FL3;
        Ndir[1].XZ() = sigma.XZ() / FL3;
        Ndir[1].XY() = sigma.XY() / FL3;
	}
	
    #ifdef LOG4CXX_PLASTICITY
    {
        LoggerPtr logger(Logger::getLogger("plasticity.SandlerDimaggio"));
        std::stringstream sout;
        sout << "<< TPZYCSandlerDimaggio::N *** \n sigma = \n" << sigma
			 << "\nI1 = " << I1 
			 << "\nJ2 = " << J2
			 << "\nSQRTJ2 = " << SQRTJ2
			 << "\nNdir = \n" << Ndir;
        LOGPZ_DEBUG(logger,sout.str().c_str());
    }
    #endif
	
	return;
}

template <class T> 
inline void TPZYCSandlerDimaggio::H(const TPZTensor<T> & sigma,const T & A, TPZVec<T> & h, int checkForcedYield) const
{

 	// the termoforce A in this case is assumed to be the
	// plastic volumetric strain itself. In fact it is not,
	// but the resultant derivatives are correct for practical purposes.
	
	//The following line evaluates L, the value of the first 
	// invariant of stresses at the intersection of the
	// shear and hardening cap yield criteria / plastic potential.
	// It is first evaluated as REAL type to avoid unnecessary
	// derivatives evaluation.
	
    #ifdef LOG4CXX_PLASTICITY
        {
          LoggerPtr logger(Logger::getLogger("plasticity.SandlerDimaggio"));
          std::stringstream sout;
          sout << ">>> TPZYCSandlerDimaggio::H *** - Hardening modulus";
          LOGPZ_INFO(logger,sout.str().c_str());
        }
    #endif
	
	REAL L_REAL, X_REAL;
	ComputeX((double)shapeFAD::val(A), X_REAL);
	LInitialGuess(X_REAL, L_REAL);
	SolveL(X_REAL, L_REAL);
    
    T I1 = sigma.I1();
	
	{//f1 - Modified Drucker-Prager as shear Yield Criterium / Plastic Potential

		h[0] = exp (I1 * T(fB) ) * T(3. * fB * fC);
		// h > 0 because plastic deformation is dilatant
	}	

	{//f2 - ellipsoidal hardening/softening cap
		T FL, X, L(L_REAL);
    	ComputeX(A, X);
		SolveL(X, L); // evaluating the derivatives of L
		ComputeF(L, FL);
		T FL2 = FL * FL;
		
		h[1] = (I1 - L) / FL2 * T ( 6. / fR / fR);
		// h <=0 because plastic deformation is compactant
	}
	
	return;
	
}

template <class T>
inline void TPZYCSandlerDimaggio::SolveL(const T & X, T & L, REAL relTol) const
{
	T F, dF, res, dRes;
	
    ComputeF(L, F);
	res = F * T(fR) + X - L;
	
	int i = 0; // ensuring evaluating at least once the Newton method
	// so that the evaluation with FAD Type and converged L evaluates the
	// derivatives
	while(
		  fabs( (double)shapeFAD::val(res) ) /
		  max( (double)fabs(shapeFAD::val(X) ), .000001 ) // avoiding division by zero
		  > relTol ||
		  i<1)
	{
        i++;
		
        ComputedF(L,dF);
        dRes = dF * T(fR) - T(1.);
		
        L -= res/dRes;
																	  
	    ComputeF(L, F);
	    res = F * T(fR) + X - L;
	}
	
}

template <class T>
inline void TPZYCSandlerDimaggio::LInitialGuess(const T & X, T & L) const
{
	T FAprox;
	ComputeF(X/T(2.), FAprox);	
	L = X + FAprox * T(fR);
}

template <class T>
inline void TPZYCSandlerDimaggio::ComputeF(const T & L, T & F) const
{
	F = L * T(fB);
	F = exp(F) * T(fC);
	F = T(fA) - F;	
}

template <class T>
inline void TPZYCSandlerDimaggio::ComputedF(const T & L, T & dF) const
{
	dF = L * T(fB);
	dF = exp(dF) * T(-fC*fB);	
}

template <class T>
inline void TPZYCSandlerDimaggio::ComputeX(const T & A, T & X) const
{
	REAL ep = - 0.99 * fW;
	
	if(shapeFAD::val(A) < ep)
	{
		T dXdep = T(fW / (ep + fW) / fD);
		X = T( log( ep/fW + 1. ) / fD );
		X = X + dXdep * (A - T(ep));
    #ifdef LOG4CXX_PLASTICITY
        {
          LoggerPtr logger(Logger::getLogger("plasticity.SandlerDimaggio"));
          std::stringstream sout;
          sout << "*** TPZYCSandlerDimaggio::ComputeX *** ##### Changing hardening equation to adjusted linear - Excessive volumetric plastic strain - Please check aftwerwards if alpha = epsP.I1() to verify if results are consistent! ######";
          LOGPZ_WARN(logger,sout.str().c_str());
        }
    #endif
	}else{
		X = A / T(fW) + T(1.);
		X = log( X ) / T(fD);
	}
	
}

//////////////////CheckConv related methods/////////////////////

inline int TPZYCSandlerDimaggio::NumCases()
{
    return 4;
}

inline void TPZYCSandlerDimaggio::LoadState(TPZFMatrix &state)
{

  int i;
  for(i=0; i<6; i++) gRefTension.fData[i] = state(i,0);

#ifdef LOG4CXX_PLASTICITY
  {
  LoggerPtr logger(Logger::getLogger("plasticity.SandlerDimaggio"));

  std::stringstream sout;
  sout << ">>> LoadState *** Tension " << state;
  LOGPZ_DEBUG(logger,sout.str().c_str());
  }
#endif
}

inline void TPZYCSandlerDimaggio::ComputeTangent(TPZFMatrix &tangent, TPZVec<REAL> &, int icase)
{

  const int nVars = 6;
  typedef TFad<nVars,double> TFAD;

  int i, j;
  TPZVec< TPZTensor < REAL > > N_Dir(2);
  REAL A = -0.05; // example epsVP value to be used with checkconv
  TPZTensor < TFAD > Sigma_FAD;
  TFAD A_FAD(A);
  TPZVec< TPZTensor < TFAD > > N_Dir_FAD(2);

  switch(icase)
  {
    case 0:
	case 1:
      //Compute N
      tangent.Redim(1,nVars);
      N(gRefTension, A, N_Dir, 0);
      for(i=0; i<nVars; i++)
          tangent(0,i) = N_Dir[icase].fData[i];
      break;
    case 2:
	case 3:
      //Compute derivatives of N
      tangent.Redim(nVars,nVars);
      gRefTension.CopyTo(Sigma_FAD);
      for(i = 0;i < nVars; i++)
         Sigma_FAD.fData[i].diff(i,nVars);
      N(Sigma_FAD, A_FAD, N_Dir_FAD, 0);
      for(i = 0; i < nVars; i++)
        for(j = 0; j < nVars; j++)
          tangent(i,j) = N_Dir_FAD[icase-2].fData[i].dx(j);
    break;

  }
#ifdef LOG4CXX_PLASTICITY
  LoggerPtr logger(Logger::getLogger("plasticity.SandlerDimaggio"));
  std::stringstream sout;
  sout << ">>> ComputeTangent *** " << tangent;
  LOGPZ_DEBUG(logger,sout.str().c_str());
#endif
}

inline void TPZYCSandlerDimaggio::Residual(TPZFMatrix &res,int icase)
{

  int i;
  const int nVars = 6;
  typedef TFad<nVars,double> TFAD;

  TPZVec< REAL > PlasticPot(2);
  REAL A = -0.05; // example epsVP value to be used with checkconv
  TPZVec< TPZTensor<REAL> > N_Dir(2);
  switch(icase)
  {
    case 0:
	case 1:
      //Compute PlasticPotential
      res.Redim(1,1);
      Compute(gRefTension, A, PlasticPot, 0);
      res(0,0) = PlasticPot[icase];
    break;
    case 2:
	case 3:
      //Compute Ndir
      res.Redim(nVars,1);
      N(gRefTension, A, N_Dir,0);
      for(i = 0; i < nVars; i++)
         res(i,0) = N_Dir[icase-2].fData[i];
    break;
  }
	
#ifdef LOG4CXX_PLASTICITY
  LoggerPtr logger(Logger::getLogger("plasticity.SandlerDimaggio"));
  std::stringstream sout;
  sout << ">>> Residual *** " << res;
  LOGPZ_DEBUG(logger,sout.str().c_str());
#endif
  
}

inline void TPZYCSandlerDimaggio::CheckConv()
{

   const int nVars = 6;

   // Creating the Sandler Dimaggio obejct
   TPZYCSandlerDimaggio YCSandlerDimaggio;

   TPZYCSandlerDimaggio::McCormicRanchSand(YCSandlerDimaggio);

   REAL multipl = 1.;
	
   TPZFNMatrix<nVars> sigma(nVars,1), Range(nVars,1);
   sigma(_XX_,0) = -0.17*multipl;
   sigma(_YY_,0) = -0.13*multipl;
   sigma(_ZZ_,0) = -0.11*multipl;
   sigma(_XY_,0) = -0.7 *multipl;
   sigma(_XZ_,0) = -0.5 *multipl;
   sigma(_ZZ_,0) = -0.3 *multipl;

   Range = sigma * (1./19.);
   TPZVec< REAL > Coefs(1,1.);
   CheckConvergence(YCSandlerDimaggio, sigma, Range, Coefs);
   
}

//////////////////CheckConv related methods/////////////////////

//////////////////Internal routines verification/////////////////

inline void TPZYCSandlerDimaggio::TestSolveL()
{
   #ifdef LOG4CXX_PLASTICITY
   LoggerPtr loggerSandlerDimaggio(Logger::getLogger("plasticity.SandlerDimaggio"));
   {
     std::stringstream sout;
     sout << "<<< TPZYCSandlerDimaggio::TestSolveL ***";
     LOGPZ_INFO(loggerSandlerDimaggio,sout.str().c_str());
   }
   #endif
	
   const int oneVar = 1;
   typedef TFad<oneVar, double> TFAD_ONE;
	
   // Creating the Sandler Dimaggio obejct
   TPZYCSandlerDimaggio YCSandlerDimaggio;

   TPZYCSandlerDimaggio::McCormicRanchSand(YCSandlerDimaggio);

   // Verifying if the SolveL routines are working
   // It needs a correct evaluation of both F and dF
	
   //Supposing a plastic volumetric strain
   REAL X, L, dF, epsVP = -.05;
   TFAD_ONE L_FAD, F_FAD;
   YCSandlerDimaggio.ComputeX(epsVP, X);

   YCSandlerDimaggio.LInitialGuess(X, L);
   L_FAD = L;
   L_FAD.diff(0, 0);
	
   YCSandlerDimaggio.ComputeF(L_FAD, F_FAD);
   YCSandlerDimaggio.ComputedF((double)shapeFAD::val(L_FAD), dF);
	
   #ifdef LOG4CXX_PLASTICITY
   {
     std::stringstream sout;
     sout << "*** TPZYCSandlerDimaggio::TestSolveL ***";
     sout << "\nDerivative evaluated with FAD (dF/dL)FAD = " << F_FAD.dx(0);
     sout << "\nDerivative evaluated explicitly (dF/dL) =  " << dF;
     sout << "\nProposed L = " << L;
     LOGPZ_DEBUG(loggerSandlerDimaggio,sout.str().c_str());
   }
   #endif
  
   YCSandlerDimaggio.SolveL(X, L);
   #ifdef LOG4CXX_PLASTICITY
   {
     std::stringstream sout;
     sout << "*** TPZYCSandlerDimaggio::TestSolveL ***";
     sout << "\nSolved   L = " << L;
     LOGPZ_DEBUG(loggerSandlerDimaggio,sout.str().c_str());
   }
   #endif
		
   //REAL multipl = 1.; // testing the shear yield criterium (f1)
   REAL multipl = 10.; // testing the elipsoidal compressive cap (f2)
	
   // Checking if NDir:(1,1,1,0,0,0) equals H
   // verifying if the TFAD derivatives of N equal the Ndir vector.
   const int sixVars = 6;
   typedef TFad<sixVars, double> TFAD_SIX;

   TPZTensor<REAL> sigma;
   TPZVec<TPZTensor<REAL> > Ndir(2);
   TPZVec<REAL> h(2);
   TPZTensor<TFAD_SIX> sigma_FAD;
   TFAD_SIX epsVP_FAD(epsVP);
   TPZVec< TFAD_SIX > PlasticPot(2);
	
   sigma.XX() = -0.17*multipl;
   sigma.YY() = -0.13*multipl;
   sigma.ZZ() = -0.11*multipl;
   sigma.XY() = -0.7 *multipl;
   sigma.YZ() = -0.5 *multipl;
   sigma.XZ() = -0.3 *multipl;
	
   sigma_FAD.XX() = sigma.XX();
   sigma_FAD.XX().diff(0,sixVars);
	
   sigma_FAD.YY() = sigma.YY();
   sigma_FAD.YY().diff(3,sixVars);
	
   sigma_FAD.ZZ() = sigma.ZZ();
   sigma_FAD.ZZ().diff(5,sixVars);
	
   sigma_FAD.XY() = sigma.XY();
   sigma_FAD.XY().diff(1,sixVars);
	
   sigma_FAD.YZ() = sigma.YZ();
   sigma_FAD.YZ().diff(4,sixVars);
	
   sigma_FAD.XZ() = sigma.XZ();
   sigma_FAD.XZ().diff(2,sixVars);
	
   YCSandlerDimaggio.N(sigma, epsVP, Ndir, 0);
   YCSandlerDimaggio.H(sigma, epsVP, h, 0);
   YCSandlerDimaggio.Compute(sigma_FAD, epsVP_FAD, PlasticPot, 0);

   #ifdef LOG4CXX_PLASTICITY
   {
	 std::stringstream sout;
     sout << "*** TPZYCSandlerDimaggio::TestSolveL ***";
     sout << "\nVerifying if H equals depsVP/dSigma = Ndir:(1 1 1 0 0 0)" ;
	
     for(int i = 0; i < NYield; i++)
     {

	    sout << "\n" << i << "-th Yield Criterium";
        sout << "\nepsVP calculated from Ndir = " << Ndir[i].XX()+Ndir[i].YY()+Ndir[i].ZZ();
        sout << "\nepsVP calculated from H =    " << h[i];
        sout << "\nVerifying if dPlasticPot/dSigma equals Ndir" ;
	    sout << "\nNdir evaluated explicitly with function N():" << Ndir[i];
        sout << "\nN evaluated through FAD evaluations: dPlasticPot/dsigma:" << PlasticPot[i];

	 }

	 LOGPZ_DEBUG(loggerSandlerDimaggio,sout.str().c_str()); 
   }
   #endif
	
}
//////////////////Internal routines verification/////////////////


inline void TPZYCSandlerDimaggio::McCormicRanchSand(TPZYCSandlerDimaggio & material)
{
   #ifdef LOG4CXX_PLASTICITY
   LoggerPtr loggerSandlerDimaggio(Logger::getLogger("plasticity.SandlerDimaggio"));
   {
      std::stringstream sout;
      sout << ">>> TPZYCSandlerDimaggio::McCormicRanchSand ***";
      LOGPZ_INFO(loggerSandlerDimaggio,sout.str().c_str());
   }
   #endif
   // setup with data from McCormic Ranch Sand
   // Dimaggio, Frank L. Sandler, Ivan S. Material model for granular soils
   // J. of the Eng. Mech. Div. vol. 97 n0 EM3 
   // pp 935-949,1971
   // OBS: stresses in ksi
   REAL A = 0.25,
        B = 0.67,
        C = 0.18,
        D = 0.67,
        R = 2.5,
        W = 0.066;
	
   material.SetUp(A, B, C, D, R, W);
}

#endif //TPZYCSANDLERDIMAGGIO_H
