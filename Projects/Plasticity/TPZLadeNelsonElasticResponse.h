// $Id: TPZLadeNelsonElasticResponse.h,v 1.17 2010-06-11 22:12:14 diogo Exp $

#ifndef TPZLADENELSONELASTICRESPONSE_H
#define TPZLADENELSONELASTICRESPONSE_H

#include "TPZTensor.h"
#include "pzreal.h"
#include "pzfmatrix.h"
#include "pzstepsolver.h"
#include "pzvec_extras.h"
#include "pzdiffmatrix.h"


#ifndef CHECKCONV
#define CHECKCONV
#include "checkconv.h"
#endif

using namespace std;

#include "tfad.h"
#include "fadType.h"
#include "pzlog.h"

#ifdef LOG4CXX_PLASTICITY
static LoggerPtr loggerPlasticity(Logger::getLogger("plasticity.plasticstep"));
#endif

/**
Calcula a tensao em funcao de deformacao (elastica)
*/
class TPZLadeNelsonElasticResponse 
{

public:

    TPZLadeNelsonElasticResponse() : fLambda(0.), fM(0.), fPoisson(0.), fPa(0.)
    { }
	
    TPZLadeNelsonElasticResponse(const TPZLadeNelsonElasticResponse & source)
    {
	   fLambda	= source.fLambda;
       fM		= source.fM;
	   fPoisson	= source.fPoisson;
	   fPa		= source.fPa;
    }

    TPZLadeNelsonElasticResponse & operator=(const TPZLadeNelsonElasticResponse & source)
    {
	   fLambda	= source.fLambda;
       fM		= source.fM;
	   fPoisson	= source.fPoisson;
	   fPa		= source.fPa;
		return *this;
    }
	
	const char * Name() const
    {
	   return "TPZLadeNelsonElasticResponse";	
    }
	
    void Print(std::ostream & out) const
    {
		out << "\n" << this->Name();
		out << "\n fLambda  = " << fLambda;
		out << "\n fM       = " << fM;
		out << "\n fPoisson = " << fPoisson;
		out << "\n fPa      = " << fPa;
    }

	
    /**
    Construtor da classe em funcao das variaveis de lame
    @param [in] Lambda Exponent in the Lade/Nelson elasticity law
    @param [in] M Multiplier constant in the Lade/Nelson elasticity law
    @param [in] poisson Poisson's coeficient
    @param [in] pa atmospheric pressure
    */
    void SetUp(REAL Lambda, REAL M, REAL poisson, REAL pa)
    {
    	fM       = M;
    	fLambda  = Lambda;
	    fPoisson = poisson;
        fPa      = pa;
    }
   
    /**
     * Computes the stress tensor based on the strain tensor
     * @param epsilon_T [in] strain tensor
     * @param sigma_T [in/out] stress tensor. It is used as the initial guess and also as output parameter
     */
    template <class T>
    void Compute(const TPZTensor<T> & epsilon_T, TPZTensor<T> & sigma_T) const;
	
    /**
     * Computes the stress tensor based on the strain tensor.
	 * REAL type specialization
     * @param epsilon_T [in] strain tensor
     * @param sigma_T [in/out] stress tensor. It is used as the initial guess and also as output parameter
     */
    void Compute(const TPZTensor<REAL> & epsilon, TPZTensor<REAL> & sigma) const;

protected:

    /**
     * Computes the stress tensor based on the strain tensor
	 * Uses a Newton's method because the YoungModulus is explicit on sigma and not epsilon
     * @param epsilon_T [in] strain tensor
     * @param sigma_T [in/out] stress tensor. It is used as the initial guess and also as output parameter
     */
    template <class T>
    void SolveSigma(const TPZTensor<T> & epsilon_T, TPZTensor<T> & sigma_T) const;
		
    template <class T, class TBASE>
    void ComputeYoung(const TPZTensor<T> & sigma, T & Young) const;

    template <class T, class TBASE>
    void ApplyElasticTensor(const T & Young, const TPZTensor<T> & Epsilon, TPZTensor<T> & sigma) const;

	template <class T, class VECTOR, class MATRIX>
    inline void ExtractTangent(const TPZTensor<T> & Res_T, VECTOR & ResVal, REAL & resnorm, MATRIX & tangent) const;
	
public:
	
    REAL fLambda;
    REAL fM;
    REAL fPoisson;
    REAL fPa;

public:
//////////////////CheckConv related methods/////////////////////

    /**
    number of types of residuals
    */
    int NumCases()
    {
      return 2;
    }

    static TPZTensor<REAL> gRefDeform;
    /**
    LoadState will keep a given state as static variable of the class
    */
    void LoadState(TPZFMatrix &state)
    {
    #ifdef LOG4CXX_PLASTICITY
        LoggerPtr logger(Logger::getLogger("plasticity.erladenelson"));
    #endif
      int i;
      for(i=0; i<6; i++) gRefDeform.fData[i] = state(i,0);
	#ifdef LOG4CXX_PLASTICITY
      std::stringstream sout;
      sout << "State " << state;
      LOGPZ_DEBUG(logger,sout.str().c_str());
	#endif
    }

    void ComputeTangent(TPZFMatrix &tangent, TPZVec<REAL> &, int icase)
    {
    #ifdef LOG4CXX_PLASTICITY
        LoggerPtr logger(Logger::getLogger("plasticity.erladenelson"));
    #endif

      const int nVars = 6;
      typedef TFad<nVars,REAL> TFAD;

      int i, j;
      TPZTensor<TFAD> epsilon_FAD, sigma_FAD;
      TFAD young_FAD;

      switch(icase)
      {
        case 1:
          //ComputeYoung
          gRefDeform.CopyTo(sigma_FAD);
          for(i=0;i<nVars;i++) 
            sigma_FAD.fData[i] . diff(i,nVars);
          tangent.Redim(1,nVars);
          ComputeYoung<TFAD, REAL>(sigma_FAD, young_FAD);
          for(i=0; i<nVars; i++)
              tangent(0,i) = young_FAD.dx(i);
          break;

        case 0:
          gRefDeform.CopyTo(epsilon_FAD);
          for(i=0;i<nVars;i++) 
            epsilon_FAD.fData[i] . diff(i,nVars);
          epsilon_FAD.Multiply(-2.,1.);// multiplying by -2 to analise the effect of previous evaluated functions
          Compute(epsilon_FAD, sigma_FAD);
          tangent.Redim(nVars,nVars);
          for(i=0; i<nVars; i++)
            for(j=0; j<nVars; j++)
              tangent(i,j) = sigma_FAD.fData[i].dx(j);
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
        LoggerPtr logger(Logger::getLogger("plasticity.erladenelson"));
    #endif
      int i;
      const int nVars = 6;
      TPZTensor< REAL > epsilon, sigma;
      REAL young;
      switch(icase)
      {
        case 1:
          //ComputeYoung
          gRefDeform.CopyTo(sigma);
          res.Redim(1,1);
          ComputeYoung<REAL , REAL>(sigma, young);
          res(0,0) = young;
        break;

        case 0:
          //
          gRefDeform.CopyTo(epsilon);
          epsilon.Multiply(-2.,1.); // multiplying by -2 to analise the effect of previous evaluated functions
          Compute(epsilon, sigma);
          res.Redim(nVars,1);
          for(i=0; i<nVars; i++)
              res(i,0) = sigma.fData[i];
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
   const int nVars = 6;

   // Creating the LadeNelsonElasticResponse obejct
   TPZLadeNelsonElasticResponse ERLadeNelson;
   // setup with data from Sacramento River Sand
   // Lade, Poul V., Kim, Moon K. Single Hardening Constitutive Model
   // for soil, Rock and Concrete. Int. J. Solids Structures, vol32 
   // No14 pp 1963-1978,1995
   ERLadeNelson.SetUp(0.28 /*Lambda*/, 900. /*M*/, 0.2 /*poisson*/, 14.7);
   TPZFNMatrix<nVars> Epsilon(nVars,1), Range(nVars,1);
   Epsilon(_XX_,0) = 0.17*0.10;
   Epsilon(_YY_,0) = 0.13*0.10;
   Epsilon(_ZZ_,0) = 0.11*0.10;
   Epsilon(_XY_,0) = 0.7 *0.10;
   Epsilon(_XZ_,0) = 0.5 *0.10;
   Epsilon(_YZ_,0) = 0.3 *0.10;

   Range = Epsilon * (1./19.);
   TPZVec< REAL > Coefs(1,1.);
   CheckConvergence(ERLadeNelson, Epsilon, Range, Coefs);
}


//////////////////CheckConv related methods/////////////////////


};


template <class T>
inline void TPZLadeNelsonElasticResponse::
        SolveSigma(const TPZTensor<T> & epsilon_T, TPZTensor<T> & sigma_T) const
{
	const int nVars = 6;
	typedef TFad< nVars, T > TFAD_SIX;
	
	TFAD_SIX Young_FAD(T (0.));
	
	TPZTensor< TFAD_SIX > sigma_FAD(Young_FAD), // initializing with zero values
	                      EEpsilon_FAD(Young_FAD),// Explicit initialization is
	                      Res_FAD(Young_FAD),// necessary because 3rd derivatives
	                      epsilon_FAD(Young_FAD);// may arise in the code
	
	TPZDiffMatrix<T> DResDSigma_T(nVars,nVars), Res_T(nVars,1), Sol_T(nVars, 1);
	
	int i;
	REAL residual, refResidual = 0.;
	EStatus status;
	
	for(i = 0 ; i < nVars; i++)
	{
		epsilon_FAD.fData[i].val() = epsilon_T.fData[i];
		sigma_FAD  .fData[i].val() = sigma_T  .fData[i];
		sigma_FAD  .fData[i].diff(i,nVars);
		//refResidual += pow(shapeFAD::val(sigma_T.fData[i]),2.);
	}
	//refResidual = max(sqrt(refResidual),1.e-6); // avoiding division by zero
	
	ComputeYoung<TFAD_SIX, T>(sigma_FAD, Young_FAD);
	ApplyElasticTensor<TFAD_SIX, T>(Young_FAD, epsilon_FAD, EEpsilon_FAD);
	//for(i = 0; i < nVars; i++)
	Res_FAD =  sigma_FAD;
	Res_FAD.Add(EEpsilon_FAD, T(-1.));
	ExtractTangent(Res_FAD, Res_T, residual, DResDSigma_T);
	
	do{
		Sol_T = Res_T;
		status = DResDSigma_T.Decompose_LU();
		status = DResDSigma_T.Substitution(& Sol_T);
		
		for(i = 0; i < nVars; i++)sigma_FAD.fData[i].val() -= Sol_T(i,0);
		
		ComputeYoung<TFAD_SIX,T>(sigma_FAD, Young_FAD);
        ApplyElasticTensor<TFAD_SIX,T>(Young_FAD, epsilon_FAD, EEpsilon_FAD);
		
		refResidual = 0.;
		for(i = 0; i < nVars; i++)refResidual += pow(shapeFAD::val(EEpsilon_FAD.fData[i]),2.);
		refResidual = sqrt(refResidual);
		
        Res_FAD =  sigma_FAD;
        Res_FAD.Add(EEpsilon_FAD, T(-1.));
		
        ExtractTangent(Res_FAD, Res_T, residual, DResDSigma_T);
	}while(residual > 1.e-6 || residual/refResidual > 1.e-12); 
	//residual > 1.e-6 ensures that the residual will be small
	// and residual/refResidual > 1.e-12 ensures that the residual is 9 orders
	// of magnitude smaller than the stress state when it converges.
	
	for(i = 0; i < nVars; i++ )
		sigma_T.fData[i] = sigma_FAD.fData[i].val();
}

inline void TPZLadeNelsonElasticResponse::
        Compute(const TPZTensor<REAL> & epsilon, TPZTensor<REAL> & sigma) const
{
	REAL YoungGuess = 1.e10;
	
	ApplyElasticTensor<REAL, REAL>(YoungGuess, epsilon, sigma); // sigma initial guess coherent with the Young guess
	
	SolveSigma(epsilon, sigma); // sigma implicitly carries the young initial guess
	// since the Young(sigma) will lead to the YoungGuess proposed value. 
	
	return;
}

template <class T>
inline void TPZLadeNelsonElasticResponse::
        Compute(const TPZTensor<T> & epsilon_T, TPZTensor<T> & sigma_T) const
{
	TPZTensor<REAL> epsilon, sigma;
	
	epsilon_T.CopyTo(epsilon);
	
	Compute(epsilon, sigma); // Solving for the correct sigma
	
    sigma.CopyTo(sigma_T);
	
	SolveSigma(epsilon_T, sigma_T); // solving again using last solution
				//as initial guess in order only to evaluate the derivatives
	
	return;
}

template <class T, class TBASE>
inline void TPZLadeNelsonElasticResponse::
         ComputeYoung(const TPZTensor<T> & sigma, T & Young) const
{
      TBASE TBase;
      REAL R = 6. * (1. + fPoisson) / (1. - 2. * fPoisson);
      T I1 = sigma.I1();
      T J2 = sigma.J2_T(TBase);
      REAL pa2 = fPa * fPa;
      //T Base = I1 * I1 / pa2 + J2 / pa2 * R;
      T Base = (I1 * I1 + J2 * TBASE(R) ) / TBASE(pa2);
      if((fLambda < 1E-12) )
      {
         Young = TBASE( fM * fPa);
         return;
      }
      REAL BaseParameter = 1.E-12;
      if((fabs(shapeFAD::val(Base)) < BaseParameter) )
      {
         // In the case the the averged stresses leads to less than
         // Base * fPa the Elastic modulus is assumed constant and
         // equal to the one at averaged stress = Base * fPa.
         // This tricky change is necessary because the Young
         // modulus increases with increasing stress state and the
         // values for very small stress states near zero, leading
         // to numerical instability (very large stress forecasts
         // in the Newton scheme at TPZPlasticStep::Sigma.
         Base = TBASE(BaseParameter);
#ifdef LOG4CXX_PLASTICITY
        {
            std::stringstream sout;
            sout << "*** ComputeYoung *** Too small Elastic Modulus - Constraining it.";
            LOGPZ_INFO(loggerPlasticity,sout.str().c_str());
        }
#endif
      }

      //Young = T( fM * fPa ) * pow( Base, fLambda );
	  Young = TBASE( fM * fPa ) * exp( TBASE(fLambda) * log( Base ) );
}

template <class T, class TBASE>
inline void TPZLadeNelsonElasticResponse::
        ApplyElasticTensor(const T & Young, const TPZTensor<T> & Epsilon, TPZTensor<T> & sigma) const
{
      // Evaluating the Lamè coefficients
	  TBASE TBase;
      T Lambda = Young * TBASE(fPoisson/((1.+fPoisson)*(1.-2.*fPoisson)));
      T Mu = Young / TBASE(2.*(1.+fPoisson));
      T trace = Epsilon.I1();
      sigma.Identity_T(TBase);
      sigma.Multiply(trace,Lambda);
      sigma.Add(Epsilon,Mu * TBASE(2.));
}

template <class T, class VECTOR, class MATRIX>
inline void TPZLadeNelsonElasticResponse::
        ExtractTangent(const TPZTensor<T> & Res_T, VECTOR & ResVal, REAL & resnorm, MATRIX & tangent) const
{
  const int nVars = 6;
  int i,j;
  typedef TFad<nVars, double> TFAD;
  // extract the values of the residual vector
  
  resnorm = 0.;
	
  for(i=0; i<nVars; i++)
  {
    ResVal(i,0) = Res_T.fData[i].val();
	resnorm += pow(shapeFAD::val(ResVal(i,0)) , 2.);
  }

  resnorm = sqrt(resnorm);

  // extract the partial derivatives to fill the tangent matrix
  tangent.Reset();
			
  for(i=0; i<nVars; i++) 
     for(j=0; j<nVars; j++) 
        tangent(i,j) = Res_T.fData[i].dx(j);

}

#endif //TPZLADEKIMELASTICRESPONSE_H

