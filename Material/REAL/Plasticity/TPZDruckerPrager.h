/**
 * @file
 */

#ifndef TZPDRUCKERPRAGER_H
#define TZPDRUCKERPRAGER_H

#include "pzlog.h"
#include "TPZPlasticStep.h"
#include "TPZYCDruckerPrager.h"
#include "TPZThermoForceA.h"
#include "TPZElasticResponse.h"
#include "pzvec_extras.h"
#include "TPZPlasticStepID.h"

#ifdef PZ_LOG
static TPZLogger loggerDrucker("plasticity.Drucker");
#endif

#define DRUCKERPARENT TPZPlasticStep<TPZYCDruckerPrager, TPZThermoForceA, TPZElasticResponse>

class TPZDruckerPrager : public DRUCKERPARENT {
public:

    enum {
        NYield = TPZYCDruckerPrager::NYield
    };

public:

    TPZDruckerPrager() : DRUCKERPARENT(), faPa(0.), fInitialEps() {
        fMaterialTensionSign = 1; // internally in this material tension is negative
        fInterfaceTensionSign = 1; // by default
    }

    TPZDruckerPrager(const TPZDruckerPrager & source) : DRUCKERPARENT(source) {
        faPa = source.faPa;
        fInitialEps = source.fInitialEps;
    }

    TPZDruckerPrager & operator=(const TPZDruckerPrager & source) {
        DRUCKERPARENT::operator=(source);
        faPa = source.faPa;
        fInitialEps = source.fInitialEps;

        return *this;
    }

    virtual const char * Name() const override {
        return "TPZDruckerPrager";
    }

    /**
         SetUp feeds all the parameters necessary to the method, distributing its values
         inside the aggregation hierarchy and computing the correct initial plasticity 
         parameter to ensure the irreversibility effect of plastic deformations.
         Elastic Mudulus:    poisson, M,    lambda
         Failure Criterion:  a,       m,    neta1
         Plastic Potential:  ksi2,    mu
         Hardening Function: C,       p
         Yield Function:     h,       alpha
         Atmospheric pression pa - to dimensionalize/adim. the stresses
     */
    void SetUp(REAL young, REAL poisson, REAL fangle, REAL coesion, REAL hardeningModulus, int InnerOuter = 0/*Inner*/) {

        DRUCKERPARENT::fYC.SetUp(fangle / 180. * M_PI, InnerOuter);
        DRUCKERPARENT::fTFA.SetUp(coesion, hardeningModulus);
        DRUCKERPARENT::fER.SetEngineeringData(young, poisson);
        //		DRUCKERPARENT::fYC.SetUp(/*phi=20*/ 20./180. * M_PI ,/*innerMCFit*/0);
        //		DRUCKERPARENT::fTFA.SetUp(/*yield- coesao inicial correspondeno a fck igual 32 Mpa */ 9.2376, /*k Modulo de hardening da coesao equivante 1 Mpa a cada 0.1% de deformacao */1000.);
        //		DRUCKERPARENT::fER.SetUp(/*young*/ 20000., /*poisson*/ 0.2);

        //		TPZTensor<REAL> nullSigma,epsA;
        //		fInitialEps = DRUCKERPARENT::GetState();

    }

    virtual void SetUp(const TPZTensor<REAL> & epsTotal)  override {
        DRUCKERPARENT::SetUp(epsTotal);
    }

    REAL YieldRadius(TPZPlasticState<REAL> state) {

        REAL radius = sqrt(2.) * DRUCKERPARENT::fTFA.Compute(state.VolHardening());
        return radius;
    }

    virtual void Print(std::ostream & out) const override {
        out << "\n" << this->Name();
        out << "\n Base Class Data:\n";
        DRUCKERPARENT::Print(out);
        out << "\nTPZDruckerPrager internal members:";
        out << "\n a*Pa = " << faPa;
        out << "\n InitialEps = " << fInitialEps;

    }

    int ClassId() const override;

    void Write(TPZStream &buf, int withclassid) const override{
        DRUCKERPARENT::Write(buf, withclassid);

        buf.Write(&faPa, 1);
        fInitialEps.Write(buf, withclassid);
    }

    void Read(TPZStream& buf, void* context) override {
        DRUCKERPARENT::Read(buf, context);

        buf.Read(&faPa, 1);
        fInitialEps.Read(buf, context);
    }

    virtual int GetNYield() const {
        return as_integer(NYield);
    }
    
public:

    // The following static members load test data from article 
    // Lade, Paul V. Kim, Moon K. Single Hardening Constitutove Model for Soil, Rock and Concrete. 
    // Int. Journal of Solid Structures, vol.32, No14. pp 1963-1978. Elsevier Science, 1994

    // Plain Concrete MPA

    static void PlainConcreteMPa(TPZDruckerPrager & material) {
        REAL poisson = 0.20;
        REAL young = 20000.;
        REAL fangle = 20.;
        REAL hardeningModulus = 1.;
        REAL coesion = 9.4;
        material.fResTol = 1.e-8;


        material.SetUp(young, poisson, fangle, coesion, hardeningModulus);

    }

    static void VeryRigidMaterial(TPZDruckerPrager & material) {
        REAL poisson = 0.20;
        REAL young = 10000000000.;
        REAL fangle = 20.;
        REAL hardeningModulus = 10000000000.;
        REAL coesion = 10000.;
        material.fResTol = 1.e-8;
        material.SetUp(young, poisson, fangle, coesion, hardeningModulus);

    }

    static void ConventionalConcrete(TPZDruckerPrager & material, int InnerOuter) {
        REAL pi = M_PI;
        REAL cohesion = 11.2033; //yield- coesao inicial correspondeno a fck igual 32 Mpa
        REAL phi = 20. / 180. * pi; //phi=20
        REAL hardening = 1000.; //Modulo de hardening da coesao equivante 1 Mpa a cada 0.1% de deformacao
        REAL young = 20000.;
        REAL poisson = 0.2;
        material.fYC.SetUp(phi, InnerOuter);
        material.fTFA.SetUp(cohesion, hardening);
        material.fER.SetEngineeringData(young, poisson);
    }

    static void TaludeMaterial(TPZDruckerPrager & material, int InnerOuter) {
        REAL pi = M_PI;
        REAL cohesion = 50.; //yield- coesao inicialem KPa
        REAL phi = 20. / 180. * pi; //phi=20
        REAL hardening = 10000.; //Modulo de hardening da coesao equivante 0.01 Kpa a cada 0.1% de deformacao
        REAL young = 20000.; //E em KPa
        REAL poisson = 0.49;
        material.fYC.SetUp(phi, InnerOuter);
        material.fTFA.SetUp(cohesion, hardening);
        material.fER.SetEngineeringData(young, poisson);
    }

private:


    /**
     * variable related to the stress translation to simulate material cohesion
     * Already contains dimensional factor
     */
    REAL faPa;

    /**
     * variable to store the plastic state related to the unstressed state
     * in a framework of a cohesion material.
     * Since LadeKim models the cohesion material with the stress axis translated
     * from a value of a, it is important to keep the related elastic deformation
     * in this variable to subtract it from the epstotal answers. It is done
     * because, intuitively, the initially nonstressed state should relate to an
     * initially undeformed state.
     */
    TPZPlasticState<REAL> fInitialEps;
};

#endif //TPZDruckerPrager_H
