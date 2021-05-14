/**
 * @file
 */

#ifndef TPZSANDLERDIMAGGIO_H
#define TPZSANDLERDIMAGGIO_H

#include "pzlog.h"
#include "TPZPlasticStep.h"
#include "TPZYCSandlerDimaggioL.h"
#include "TPZYCSandlerDimaggioL2.h"
#include "TPZSandlerDimaggioThermoForceA.h"
#include "TPZElasticResponse.h"
#include "pzvec_extras.h"
#include "TPZPlasticStepID.h"
#include "TPZPlasticStepTranslator.h"
#include "TPZYCSandlerDimaggioLTranslator.h"
#include "TPZSandlerDimaggioThermoForceATranslator.h"
#include "TPZElasticResponseTranslator.h"
#include "TPZYCSandlerDimaggioL2Translator.h"

#define SANDLERDIMAGGIOSTEP1 TPZPlasticStep<TPZYCSandlerDimaggioL, TPZSandlerDimaggioThermoForceA, TPZElasticResponse>
#define SANDLERDIMAGGIOSTEP1TRANSLATOR TPZPlasticStepTranslator<TPZYCSandlerDimaggioLTranslator, TPZSandlerDimaggioThermoForceATranslator, TPZElasticResponseTranslator>
#define SANDLERDIMAGGIOSTEP2 TPZPlasticStep<TPZYCSandlerDimaggioL2, TPZSandlerDimaggioThermoForceA, TPZElasticResponse>
#define SANDLERDIMAGGIOSTEP2TRANSLATOR TPZPlasticStepTranslator<TPZYCSandlerDimaggioL2Translator, TPZSandlerDimaggioThermoForceATranslator, TPZElasticResponseTranslator>


template<class SANDLERDIMAGGIOPARENT>
class TPZSandlerDimaggio : public SANDLERDIMAGGIOPARENT  {

public:

  enum {NYield = TPZYCSandlerDimaggio::NYield};

public:

    TPZSandlerDimaggio(REAL alpha=0./*-1.e-9*/):SANDLERDIMAGGIOPARENT(alpha) // avoiding nan
    {
		this->fMaterialTensionSign = 1;//Quando este numero for positivo carremento de compressao e negativo
        
    }

    TPZSandlerDimaggio(const TPZSandlerDimaggio & source):SANDLERDIMAGGIOPARENT(source)
    {
		this->fMaterialTensionSign = 1;
    }

    TPZSandlerDimaggio & operator=(const TPZSandlerDimaggio & source)
    {
		this->fMaterialTensionSign = 1;
		SANDLERDIMAGGIOPARENT::operator=(source);
		
		return *this;
    }
	
	virtual const char * Name() const override
	{
	   return "TPZSandlerDimaggio";	
	}
	
	virtual void Print(std::ostream & out) const override
	{
		out << "\n" << this->Name();
		out << "\n Base Class Data:\n";
		SANDLERDIMAGGIOPARENT::Print(out);
		out << "\nTPZSandlerDimaggio internal members: None";		
	}
	
    int ClassId() const override;

    
    void Write(TPZStream &buf, int withclassid) const override{
	   SANDLERDIMAGGIOPARENT::Write(buf, withclassid);
	   // fPlasticMem does not need to be stored
			
	}

    void Read(TPZStream& buf, void* context) override {
	   SANDLERDIMAGGIOPARENT::Read(buf, context);
	   this->fPlasticMem.Resize(0);
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
        SANDLERDIMAGGIOPARENT::fN.m_hardening = this->fYC.InitialDamage();        
        SANDLERDIMAGGIOPARENT::fER.SetEngineeringData(E, poisson);
    }
    virtual void SetUp(const TPZTensor<REAL> & epsTotal)  override {
        SANDLERDIMAGGIOPARENT::SetUp(epsTotal);
    }

    /**
    Retrieve the plastic state variables
    */
    virtual TPZPlasticState<REAL> GetState() const override
    {
        return SANDLERDIMAGGIOPARENT::GetState();
    }

    /**
    Computes the strain tensor as a function of the stress state.
    This function returns the inverse of function void Sigma(...) using a Newton's scheme.
    @param[in] sigma stress tensor
    @param[out] epsTotal deformation tensor
    */
    virtual void ApplyLoad(const TPZTensor<REAL> & sigma, TPZTensor<REAL> &epsTotal) override
    {
       SANDLERDIMAGGIOPARENT::ApplyLoad_Internal(sigma, epsTotal);
    }

    /**
    * Load the converged solution, updating the damage variables
    */
    virtual void ApplyStrain(const TPZTensor<REAL> &epsTotal) override
    {
        SANDLERDIMAGGIOPARENT::ApplyStrain_Internal(epsTotal);
    }

	/**
    * Imposes the specified strain tensor and performs plastic integration when necessary.
	*
    */
    virtual void ApplyStrainComputeDep(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma, TPZFMatrix<REAL> &Dep) override
	{

		SANDLERDIMAGGIOPARENT::ApplyStrainComputeDep_Internal(epsTotal, sigma, Dep);
		
	}
	
    virtual void ApplyStrainComputeSigma(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma, TPZFMatrix<REAL> * tangent = NULL) override
	{
        bool require_tangent_Q = true;
        if (!tangent) {
            require_tangent_Q = false;
        }
        
#ifdef PZDEBUG
        // Check for required dimensions of tangent
        if (!(tangent->Rows() == 6 && tangent->Cols() == 6)) {
            std::cerr << "Unable to compute the tangent operator. Required tangent array dimensions are 6x6." << std::endl;
            DebugStop();
        }
#endif
        
		SANDLERDIMAGGIOPARENT::ApplyStrainComputeSigma_Internal(epsTotal, sigma);
	}

    /**
    return the value of the yield functions for the given deformation
     * @param[in] epsTotal deformation tensor (total deformation
     * @param[out] phi vector of yield functions
    */
    virtual void Phi(const TPZTensor<REAL> &epsTotal, TPZVec<REAL> &phi) const override
    {
        SANDLERDIMAGGIOPARENT::Phi_Internal(epsTotal, phi);
    }
    
    virtual int GetNYield() const {
        return as_integer(NYield);
    }
	
protected:
	/**
	* Evaluates the constitutive matrix (DSigma/DEpsT) based on the data from the plastic
	* integration history without modifying it.
	*
	* @param [out] sigma resultant stress tensor
	* @param [out] Dep Incremental constitutive relation
    */
    virtual void ComputeDep(TPZTensor<REAL> & sigma, TPZFMatrix<REAL> &Dep) override
	{
		const int nyield = this->fYC.NYield;
		
		SANDLERDIMAGGIOPARENT::ComputeDep(sigma, Dep);

		TPZManVector<REAL,3> plastifLen(nyield, 0.);
		int n = this->fPlasticMem.NElements();
		REAL deltaAlpha = fabs(this->fPlasticMem[n-1].m_elastoplastic_state.m_hardening - this->fPlasticMem[1].m_elastoplastic_state.m_hardening); 
		
		this->IntegrationOverview(plastifLen);
		
		// if the plastification ocurred mainly in the first (0th) yield function, which is
		// perfectly plastic, then the stiffness matrix may be singular.
		if(  (plastifLen[0] > 0.9 && plastifLen[1] < 0.1) ||
		     (plastifLen[0] < 1.e-10 && plastifLen[1] > 0.9 && deltaAlpha < 1.e-10) )
		{
			TPZFNMatrix<6*6> D(6,6,0.);
			
			SANDLERDIMAGGIOPARENT::fER.De(D);
			
			Dep.ZAXPY(0.01, D);
			
			#ifdef PZ_LOG
				{
				  TPZLogger logger("plasticity.SandlerDimaggio");
				  std::stringstream sout;
				  sout << "*** TPZYCSandlerDimaggio::ComputeDep *** Superimposing a fraction of the Elastic Stiffness Matrix on a perfectly plastic load";
				  std::cout << "\nfPlasticLen = " << plastifLen << " deltaAlpha = " << deltaAlpha;
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
       #ifdef PZ_LOG
       TPZLogger loggerSandlerDimaggio("plasticity.SandlerDimaggio");
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
       #ifdef PZ_LOG
       TPZLogger loggerSandlerDimaggio("plasticity.SandlerDimaggio");
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
       #ifdef PZ_LOG
       TPZLogger loggerSandlerDimaggio("plasticity.SandlerDimaggio");
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
       #ifdef PZ_LOG
       TPZLogger loggerSandlerDimaggio("plasticity.SandlerDimaggio");
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
	
    static void UncDeepSandTest(TPZSandlerDimaggio & material)
    {
#ifdef PZ_LOG
        TPZLogger loggerSandlerDimaggio("plasticity.SandlerDimaggio");
        {
            std::stringstream sout;
            sout << ">>> TPZSandlerDimaggio::Unconsolidated Deep Sandstone Reservoir ***";
            LOGPZ_INFO(loggerSandlerDimaggio,sout.str().c_str());
        }
#endif
		
        REAL E = 1305, //ksi (9000MPa)
        poisson = 0.25;
        
//        material.fER.SetUp(E, poisson);
        
        REAL A = 10.,//2.61, //ksi  = 18 in MPa
        B = 0.169, // ksi   0.0245 in MPa
        C = 2.57, // ksi   = 17.7  in MPa
        D = 0.05069, //  = 0.00735 in MPa
        R = 1.5,
        W = 0.0908;
        
//        material.fYC.SetUp(A, B, C, D, R, W);
        material.SetUp(poisson, E, A, B, C, R, D, W);
	}
	
	static void UncDeepSandResPSI(TPZSandlerDimaggio & material)
    {
       #ifdef PZ_LOG
       TPZLogger loggerSandlerDimaggio("plasticity.SandlerDimaggio");
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
       #ifdef PZ_LOG
       TPZLogger loggerSandlerDimaggio("plasticity.SandlerDimaggio");
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
       #ifdef PZ_LOG
       TPZLogger loggerSandlerDimaggio("plasticity.SandlerDimaggio");
       {
          std::stringstream sout;
          sout << ">>> TPZSandlerDimaggio::PRSMat MPa ***";
          LOGPZ_INFO(loggerSandlerDimaggio,sout.str().c_str());
	   }
       #endif
		
	   REAL E = 29269.,
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

template<class SANDLERDIMAGGIOPARENT>
int TPZSandlerDimaggio<SANDLERDIMAGGIOPARENT>::ClassId() const{
    return Hash("TPZSandlerDimaggio") ^ SANDLERDIMAGGIOPARENT::ClassId() << 1;
}

#endif //TPZSANDLERDIMAGGIO_H
