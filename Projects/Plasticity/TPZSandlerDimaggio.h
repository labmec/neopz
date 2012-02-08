// $Id: TPZSandlerDimaggio.h,v 1.19 2009-12-10 23:05:00 erick Exp $

#ifndef TPZSANDLERDIMAGGIO_H
#define TPZSANDLERDIMAGGIO_H

#include "pzlog.h"
#include "TPZPlasticStep.h"
#include "TPZYCSandlerDimaggio.h"
#include "TPZSandlerDimaggioThermoForceA.h"
#include "TPZElasticResponse.h"
#include "pzvec_extras.h"
#include "pzsave.h"
#include "TPZPlasticStepID.h"

#define SANDLERDIMAGGIOPARENT TPZPlasticStep<TPZYCSandlerDimaggio, TPZSandlerDimaggioThermoForceA, TPZElasticResponse>

class TPZSandlerDimaggio : public SANDLERDIMAGGIOPARENT, public TPZSaveable  {

public:

  enum {NYield = TPZYCSandlerDimaggio::NYield};

public:

    TPZSandlerDimaggio(REAL alpha=0./*-1.e-9*/):SANDLERDIMAGGIOPARENT(alpha) // avoiding nan
    {
		fMaterialTensionSign = 1;
    }

    TPZSandlerDimaggio(const TPZSandlerDimaggio & source):SANDLERDIMAGGIOPARENT(source)
    {
		fMaterialTensionSign = 1;
    }

    TPZSandlerDimaggio & operator=(const TPZSandlerDimaggio & source)
    {
		fMaterialTensionSign = 1;
		SANDLERDIMAGGIOPARENT::operator=(source);
		
		return *this;
    }
	
	virtual const char * Name() const
	{
	   return "TPZSandlerDimaggio";	
	}
	
	virtual void Print(std::ostream & out) const
	{
		out << "\n" << this->Name();
		out << "\n Base Class Data:\n";
		SANDLERDIMAGGIOPARENT::Print(out);
		out << "\nTPZSandlerDimaggio internal members: None";		
	}
	
	virtual int ClassId() const
	{
		return TPZSANDLERDIMAGGIO_ID;	
	}
	
	virtual void Write(TPZStream &buf, int withclassid)
	{
	   TPZSaveable::Write(buf, withclassid);
		
	   buf. Write(&fYC.fA, 1);
	   buf. Write(&fYC.fB, 1);
	   buf. Write(&fYC.fC, 1);
	   buf. Write(&fYC.fD, 1);
	   buf. Write(&fYC.fR, 1);
	   buf. Write(&fYC.fW, 1);	
		
	   buf. Write(&fER.fLambda, 1);
	   buf. Write(&fER.fMu, 1);	

	   buf. Write(&fResTol, 1);
	   buf. Write(&fIntegrTol, 1);
	   buf. Write(&fMaxNewton, 1);
	   buf. Write(&fMinLambda, 1);
		
	   buf. Write(&fN.fEpsT.fData[0], 6);
	   buf. Write(&fN.fEpsP.fData[0], 6);
	   buf. Write(&fN.fAlpha, 1);
		
	   // fPlasticMem does not need to be stored
			
	}

	virtual void Read(TPZStream &buf, void *context)
	{
	   TPZSaveable::Read(buf, context);
		
	   buf. Read(&fYC.fA, 1);
	   buf. Read(&fYC.fB, 1);
	   buf. Read(&fYC.fC, 1);
	   buf. Read(&fYC.fR, 1);
	   buf. Read(&fYC.fW, 1);	
		
	   buf. Read(&fER.fLambda, 1);
	   buf. Read(&fER.fMu, 1);	
		
	   buf. Read(&fResTol, 1);
	   buf. Read(&fIntegrTol, 1);
	   buf. Read(&fMaxNewton, 1);
	   buf. Read(&fMinLambda, 1);
		
	   buf. Read(&fN.fEpsT.fData[0], 6);
	   buf. Read(&fN.fEpsP.fData[0], 6);
	   buf. Read(&fN.fAlpha, 1);
		
	   fPlasticMem.Resize(0);
	}	

    /**
    SetUp feeds all the parameters necessary to the method, distributing its values
    inside the aggregation hierarchy and computing the correct initial plasticity 
    parameter to ensure the irreversibility effect of plastic deformations.
    Elastic Mudulus:    E, poisson
    Failure Criterion:  no parameters
    Plastic Potential:  A, B, C, R
    Hardening Function: D, W
    Yield Function:     associative
    */
    void SetUp(REAL poisson, REAL E,
               REAL A, REAL B, REAL C, REAL R,
               REAL D, REAL W)
    {
       SANDLERDIMAGGIOPARENT::fYC.SetUp(A, B, C, D, R, W);
       SANDLERDIMAGGIOPARENT::fER.SetUp(E, poisson);
    }


    /**
    Retrieve the plastic state variables
    */
    virtual const TPZPlasticState<REAL> GetState() const
    {
        return SANDLERDIMAGGIOPARENT::GetState();
    }

    /**
    Computes the strain tensor as a function of the stress state.
    This function returns the inverse of function void Sigma(...) using a Newton's scheme.
    @param [in] sigma stress tensor
    @param [out] epsTotal deformation tensor
    */
    virtual void ApplyLoad(const TPZTensor<REAL> & sigma, TPZTensor<REAL> &epsTotal)
    {
       SANDLERDIMAGGIOPARENT::ApplyLoad_Internal(sigma, epsTotal);
    }

    /**
    * Load the converged solution, updating the damage variables
    */
    virtual void ApplyStrain(const TPZTensor<REAL> &epsTotal)
    {
        SANDLERDIMAGGIOPARENT::ApplyStrain_Internal(epsTotal);
    }

	/**
    * Imposes the specified strain tensor and performs plastic integration when necessary.
	*
    */
    virtual void ApplyStrainComputeDep(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma, TPZFMatrix &Dep)
	{

		SANDLERDIMAGGIOPARENT::ApplyStrainComputeDep_Internal(epsTotal, sigma, Dep);
		
	}
	
    virtual void ApplyStrainComputeSigma(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma)
	{
		SANDLERDIMAGGIOPARENT::ApplyStrainComputeSigma_Internal(epsTotal, sigma);
	}

    /**
    return the value of the yield functions for the given deformation
     * @param [in] deform deformation tensor (total deformation
     * @param [out] phi vector of yield functions
    */
    virtual void Phi(const TPZTensor<REAL> &epsTotal, TPZVec<REAL> &phi) const
    {
        SANDLERDIMAGGIOPARENT::Phi_Internal(epsTotal, phi);
    }
	
protected:
	/**
	* Evaluates the constitutive matrix (DSigma/DEpsT) based on the data from the plastic
	* integration history without modifying it.
	*
	* @param [out] sigma resultant stress tensor
	* @param [out] Dep Incremental constitutive relation
    */
    virtual void ComputeDep(TPZTensor<REAL> & sigma, TPZFMatrix &Dep)
	{
		const int nyield = fYC.NYield;
		
		SANDLERDIMAGGIOPARENT::ComputeDep(sigma, Dep);

		TPZVec<REAL> plastifLen(nyield, 0.);
		int n = fPlasticMem.NElements();
		REAL deltaAlpha = fabs(fPlasticMem[n-1].fPlasticState.fAlpha - fPlasticMem[1].fPlasticState.fAlpha); 
		
		IntegrationOverview(plastifLen);
		
		// if the plastification ocurred maily in the first (0th) yield function, which is
		// perfectly plastic, then the stiffness matrix may be singular.
		if(  (plastifLen[0] > 0.9 && plastifLen[1] < 0.1) ||
		     (plastifLen[0] < 1.e-10 && plastifLen[1] > 0.9 && deltaAlpha < 1.e-10) )
		{
			TPZFNMatrix<6*6> D(6,6,0.);
			
			SANDLERDIMAGGIOPARENT::fER.ElasticMat(D);
			
			Dep.ZAXPY(0.01, D);
			
			#ifdef LOG4CXX_PLASTICITY
				{
				  LoggerPtr logger(Logger::getLogger("plasticity.SandlerDimaggio"));
				  std::stringstream sout;
				  sout << "*** TPZYCSandlerDimaggio::ComputeDep *** Superimposing a fraction of the Elastic Stiffness Matrix on a perfectly plastic load";
				  cout << "\nfPlasticLen = " << plastifLen << " deltaAlpha = " << deltaAlpha;
				  LOGPZ_WARN(logger,sout.str().c_str());
				}
			#endif
		}
	
	}

public:


// The following static members load test data from article 
// setup with data from McCormic Ranch Sand
// Dimaggio, Frank L. Sandler, Ivan S. Material model for granular soils
// J. of the Eng. Mech. Div. vol. 97 n0 EM3 
// pp 935-949,1971

static void McCormicRanchSand(TPZSandlerDimaggio & material)
    {
       #ifdef LOG4CXX_PLASTICITY
       LoggerPtr loggerSandlerDimaggio(Logger::getLogger("plasticity.SandlerDimaggio"));
       {
          std::stringstream sout;
          sout << ">>> TPZSandlerDimaggio::McCormicRanchSand ***";
          LOGPZ_INFO(loggerSandlerDimaggio,sout.str().c_str());
	   }
       #endif
		
	   REAL E = 100, //ksi = 689 MPa
			poisson = /*0.25;*/ 0.25; // In the 1971 article, although
			          // the authors documented a poisson coeff. of 0.25
					  // they seem to have used a poisson of 0.40, as
					  // calculated in the unloading cycle of the uniaxial
					  // strain test.
       TPZYCSandlerDimaggio::McCormicRanchSand(material.fYC);
       material.fER.SetUp(E, poisson);
		
	   material.fResTol = 1.e-12;
	   material.fIntegrTol = 1.e-6;

    }

    static void McCormicRanchSandMod(TPZSandlerDimaggio & material)
    {
       #ifdef LOG4CXX_PLASTICITY
       LoggerPtr loggerSandlerDimaggio(Logger::getLogger("plasticity.SandlerDimaggio"));
       {
          std::stringstream sout;
          sout << ">>> TPZSandlerDimaggio::McCormicRanchSand ***";
          LOGPZ_INFO(loggerSandlerDimaggio,sout.str().c_str());
	   }
       #endif
		
	   REAL E = 100, //ksi
			poisson = /*0.25;*/ 0.40; // In the 1971 article, although
			          // the authors documented a poisson coeff. of 0.25
					  // they seem to have used a poisson of 0.40, as
					  // calculated in the unloading cycle of the uniaxial
					  // strain test.
       TPZYCSandlerDimaggio::McCormicRanchSand(material.fYC);
       material.fER.SetUp(E, poisson);

    }

    static void McCormicRanchSandMod2(TPZSandlerDimaggio & material)
    {
       #ifdef LOG4CXX_PLASTICITY
       LoggerPtr loggerSandlerDimaggio(Logger::getLogger("plasticity.SandlerDimaggio"));
       {
          std::stringstream sout;
          sout << ">>> TPZSandlerDimaggio::McCormicRanchSand ***";
          LOGPZ_INFO(loggerSandlerDimaggio,sout.str().c_str());
	   }
       #endif
		
	   REAL E = 100, //ksi
			poisson = /*0.25;*/ 0.40; // In the 1971 article, although
			          // the authors documented a poisson coeff. of 0.25
					  // they seem to have used a poisson of 0.40, as
					  // calculated in the unloading cycle of the uniaxial
					  // strain test.

       material.fER.SetUp(E, poisson);

	   REAL A = 0.25,
        B = 0.67,
        C = 0.18,
        D = 0.67 / 2., //letting the material behave stronger to enable higher loads without too much volumetric plastic strain 
        R = 2.5,
        W = 0.066;
	
       material.fYC.SetUp(A, B, C, D, R, W);
	}
	
   static void UncDeepSandRes(TPZSandlerDimaggio & material)
    {
       #ifdef LOG4CXX_PLASTICITY
       LoggerPtr loggerSandlerDimaggio(Logger::getLogger("plasticity.SandlerDimaggio"));
       {
          std::stringstream sout;
          sout << ">>> TPZSandlerDimaggio::Unconsolidated Deep Sandstone Reservoir ***";
          LOGPZ_INFO(loggerSandlerDimaggio,sout.str().c_str());
	   }
       #endif
		
	   REAL E = 1305, //ksi (9000MPa)
	   poisson = 0.25;

       material.fER.SetUp(E, poisson);

        REAL A = 2.61, //ksi  = 18 in MPa
        B = 0.169, // ksi   0.0245 in MPa
        C = 2.57, // ksi   = 17.7  in MPa
        D = 0.05069, //  = 0.00735 in MPa
        R = 1.5,
        W = 0.0908;
	
       material.fYC.SetUp(A, B, C, D, R, W);
	}
	
	static void UncDeepSandResPSI(TPZSandlerDimaggio & material)
    {
       #ifdef LOG4CXX_PLASTICITY
       LoggerPtr loggerSandlerDimaggio(Logger::getLogger("plasticity.SandlerDimaggio"));
       {
          std::stringstream sout;
          sout << ">>> TPZSandlerDimaggio::Unconsolidated Deep Sandstone Reservoir ***";
          LOGPZ_INFO(loggerSandlerDimaggio,sout.str().c_str());
	   }
       #endif
		
	   REAL E = 1305000, //psi (9000MPa)
	   poisson = 0.25;

       material.fER.SetUp(E, poisson);

        REAL A = 2610, //psi  = 18 in MPa
        B = 0.000169, // psi   0.0245 in MPa
        C = 2570, // psi    = 17.7  in MPa
        D = 0.00005069, //  = 0.00735 in MPa
        R = 1.5,
        W = 0.0908;
	
       material.fYC.SetUp(A, B, C, D, R, W);
	}

	static void UncDeepSandResMPa(TPZSandlerDimaggio & material)
    {
       #ifdef LOG4CXX_PLASTICITY
       LoggerPtr loggerSandlerDimaggio(Logger::getLogger("plasticity.SandlerDimaggio"));
       {
          std::stringstream sout;
          sout << ">>> TPZSandlerDimaggio::Unconsolidated Deep Sandstone Reservoir ***";
          LOGPZ_INFO(loggerSandlerDimaggio,sout.str().c_str());
	   }
       #endif
		
	   REAL E = 9000, 
	   poisson = 0.25;

       material.fER.SetUp(E, poisson);

        REAL A = 18, 
        B = 0.0245, 
        C = 17.7, 
        D = 0.00735, 
        R = 1.5,
        W = 0.0908;
	
       material.fYC.SetUp(A, B, C, D, R, W);
	}
	
	static void PRSMatMPa(TPZSandlerDimaggio & material)
    {
       #ifdef LOG4CXX_PLASTICITY
       LoggerPtr loggerSandlerDimaggio(Logger::getLogger("plasticity.SandlerDimaggio"));
       {
          std::stringstream sout;
          sout << ">>> TPZSandlerDimaggio::PRSMat MPa ***";
          LOGPZ_INFO(loggerSandlerDimaggio,sout.str().c_str());
	   }
       #endif
		
	   REAL E = 29269, 
	   poisson = 0.203;

       material.fER.SetUp(E, poisson);

        REAL A = 116.67, 
        B = 0.0036895, 
        C = 111.48, 
        D = 0.018768, 
        R = 0.91969,
        W = 0.006605;
	
       material.fYC.SetUp(A, B, C, D, R, W);
	}

public:

};


#endif //TPZSANDLERDIMAGGIO_H
