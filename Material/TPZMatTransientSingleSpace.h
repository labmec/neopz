/**
 * @file
 * @brief Contains the TPZMatTransientSingleSpace class which implementes an Euler integrator for discontinuous materials.
 */

#ifndef TRANSIENTMATH
#define TRANSIENTMATH

#include "TPZMaterialDataT.h"

template<class T>
class TPZBndCondT;
/**
 * @ingroup matsinglespace
 * @brief Defines an implicit Euler time integrator for using with single space discontinuous materials.
 * Weak statement is supposed to be Integral[(un+1 - un)/deltaT * v, Omega] + Bilinear Form = Linear Form
 * This class implements only Integral[(un+1 - un)/deltaT * v, Omega]. Bilinear and linear form must be implemented in the material class.
 */

template<class TBASEMAT>
class TPZMatTransientSingleSpace : public TBASEMAT {
    public:

	/** @brief Class constructor */
	TPZMatTransientSingleSpace(int nummat, int dim, REAL TimeStep);

    /** @brief Sets integral scheme as an explicit Euler */
	void SetExplicit();
	
	/** √Sets integral scheme as an implicit Euler */
	void SetImplicit();
	
	void Contribute(const TPZMaterialDataT<STATE> &data,
                    REAL weight,
                    TPZFMatrix<STATE> &ek,
                    TPZFMatrix<STATE> &ef) override;
	
	void ContributeBC(const TPZMaterialDataT<STATE> &data,
                      REAL weight,
                      TPZFMatrix<STATE> &ek,
                      TPZFMatrix<STATE> &ef,
                      TPZBndCondT<STATE> &bc) override;
	
	void ContributeInterface(const TPZMaterialDataT<STATE> &data,
                             const TPZMaterialDataT<STATE> &dataleft,
                             const TPZMaterialDataT<STATE> &dataright,
                             REAL weight,
                             TPZFMatrix<STATE> &ek,
                             TPZFMatrix<STATE> &ef) override;
	
	void ContributeBCInterface(const TPZMaterialDataT<STATE> &data,
                               const TPZMaterialDataT<STATE> &dataleft,
                               REAL weight,
                               TPZFMatrix<STATE> &ek,
                               TPZFMatrix<STATE> &ef,
                               TPZBndCondT<STATE> &bc) override;
	
	/** @brief Set material to compute only Integral[- un/deltaT * v, Omega] */
	void SetLastState();
	
	/** @brief Set material to compute Integral[un+1/deltaT * v, Omega] + Bilinear Form = Linear Form  */
	void SetCurrentState();
	
	/** @brief Set material to compute ek = Integral[phi_i phi_j, Omega]/deltaT */
	void SetMassMatrix();
	
	/** @brief Set material to compute ef = Linear Form - Bilinear Form(u) = F -ku */ 
	void SetFluxOnly();
	
	/** @brief Define time step DeltaT */
	void SetTimeStep(STATE TimeStep);
	
	/** @brief Returns time step value. */
	STATE TimeStep();
	
	/** @brief Indicates if the material requires the solution to compute Contribute */
	/** 
	 * By default its value is true, but it can be set as false by derived material classes \n
	 * to increase the performance of method TPZCompEl::CalcStiff
	 */
	virtual bool NeedsSolutionToContribute(){
		return true;
	}
	
	/** @brief Indicates if the material requires the global coordinate X to compute Contribute */
	/**
	 * By default its value is true, but it can be set as false by derived material classes \n
	 * to increase the performance of method TPZCompEl::CalcStiff
	 */
	virtual bool NeedsXCoord(){
		return (this->fStep != ELast);
	}

    int ClassId() const override;
protected:
	
	enum ETemporalScheme{EImplicit = 1, EExplicit = 2};
	
	ETemporalScheme fTemporalIntegrator;
	
	enum STEPS{ENone = -1, ELast = 0, ECurrent = 1, EMassMatrix = 2, EFluxOnly = 3};
	
	STEPS fStep;
	
	STATE fTimeStep;
	
	virtual void ContributeSolutionRhs(const TPZVec<STATE> &sol,
                                       const TPZFMatrix<REAL> &phi,
                                       REAL weight, TPZFMatrix<STATE> &ef);
	
	virtual void ContributeTangent(const TPZFMatrix<REAL> &phi,
                                   REAL weight, TPZFMatrix<STATE> &ek);
};

template<class TBASEMAT>
inline void TPZMatTransientSingleSpace< TBASEMAT >::SetLastState(){
	this->fStep = ELast;
}

template<class TBASEMAT>
inline void TPZMatTransientSingleSpace< TBASEMAT >::SetCurrentState(){
	this->fStep = ECurrent;
}

template<class TBASEMAT>
inline void TPZMatTransientSingleSpace< TBASEMAT >::SetMassMatrix(){
	this->fStep = EMassMatrix;
}

template<class TBASEMAT>
inline void TPZMatTransientSingleSpace< TBASEMAT >::SetFluxOnly(){
	this->fStep = EFluxOnly;
}

template<class TBASEMAT>
inline void TPZMatTransientSingleSpace< TBASEMAT >::SetTimeStep(REAL TimeStep){
	this->fTimeStep = TimeStep;
}

template<class TBASEMAT>
inline REAL TPZMatTransientSingleSpace< TBASEMAT >::TimeStep(){
	return this->fTimeStep;
}


template<class TBASEMAT>
int TPZMatTransientSingleSpace<TBASEMAT>::ClassId() const{
    return Hash("TPZMatTransientSingleSpace") ^ TBASEMAT::ClassId() << 1;
}
template<class TBASEMAT>
void TPZMatTransientSingleSpace< TBASEMAT >::SetExplicit() {
	this->fTemporalIntegrator = EExplicit;
}
template<class TBASEMAT>
void TPZMatTransientSingleSpace< TBASEMAT >::SetImplicit() {
	this->fTemporalIntegrator = EImplicit;
}

template<class TBASEMAT>
TPZMatTransientSingleSpace< TBASEMAT >::TPZMatTransientSingleSpace(int nummat, int dim, REAL TimeStep):TPZRegisterClassId(&TPZMatTransientSingleSpace::ClassId),
TBASEMAT(nummat, dim) {
	this->SetTimeStep(TimeStep);
}

template<class TBASEMAT>
void TPZMatTransientSingleSpace< TBASEMAT >::Contribute(const TPZMaterialDataT<STATE> &data,
                                                  REAL weight,
                                                  TPZFMatrix<STATE> &ek,
                                                  TPZFMatrix<STATE> &ef) {
	
	// Mostly for implicit
    int numbersol = data.sol.size();
    if (numbersol != 1) {
        DebugStop();
    }

	if (this->fStep == ECurrent){
		TBASEMAT::Contribute(data,weight,ek,ef);
		ContributeSolutionRhs(data.sol[0], data.phi, weight, ef);
		ContributeTangent(data.phi, weight, ek);
		return;
	}
	
	if (this->fStep == ELast){
		this->ContributeSolutionRhs(data.sol[0], data.phi, weight, ef);
		return;
	}
	
	// Mostly for explicit
	if (this->fStep == EMassMatrix){
		this->ContributeTangent(data.phi, weight, ek);
		return;
	}
	
	if (this->fStep == EFluxOnly){ //Calcula ef = F-ku
		TBASEMAT::Contribute(data,weight,ek,ef);
		return;
	}
	
	
	PZError << "ERROR! " << __PRETTY_FUNCTION__ << " at LINE " << __LINE__ << std::endl;
	
}

template<class TBASEMAT>
void TPZMatTransientSingleSpace< TBASEMAT >::ContributeBC(const TPZMaterialDataT<STATE> &data,
                                                    REAL weight,
                                                    TPZFMatrix<STATE> &ek,
                                                    TPZFMatrix<STATE> &ef,
                                                    TPZBndCondT<STATE> &bc) {
	// Mostly for implicit
	if (this->fStep == ECurrent){
		TBASEMAT::ContributeBC(data,weight,ek,ef,bc);
		return;
	}
	
	if (this->fStep == ELast){
		return;
	}
	
	
	// Mostly for explicit
	if (this->fStep == EMassMatrix){
		TPZFNMatrix<1000,STATE> fakeef(ek.Rows(),1,0.);
		TBASEMAT::ContributeBC(data,weight,ek,fakeef,bc);
		return;
	}
	if (this->fStep == EFluxOnly){ //Calcula ef = F-ku
		TPZFNMatrix<1000> fakeef(ef.Rows(),ef.Rows(),0.);
		TBASEMAT::ContributeBC(data,weight,ek,ef,bc);
		return;
	}
	
	
	PZError << "ERROR! " << __PRETTY_FUNCTION__ << " at LINE " << __LINE__ << std::endl;
}

template<class TBASEMAT>
void TPZMatTransientSingleSpace< TBASEMAT >::ContributeInterface(
    const TPZMaterialDataT<STATE> &data,
    const TPZMaterialDataT<STATE> &dataleft,
    const TPZMaterialDataT<STATE> &dataright,
    REAL weight,
    TPZFMatrix<STATE> &ek,
    TPZFMatrix<STATE> &ef)
{
	
	// Mostly for implicit
	if (this->fStep == ECurrent){
		TBASEMAT::ContributeInterface(data,dataleft,dataright, weight, ek, ef);
		return;
	}
	
	if (this->fStep == ELast){
		return;
	}
	
	// Mostly for explicit
	if (this->fStep == EMassMatrix){
		return;
	}
	if (this->fStep == EFluxOnly){ //Calcula ef = F-ku
		TBASEMAT::ContributeInterface(data,dataleft,dataright, weight, ek, ef);
		return;
	}
	
	
	PZError << "ERROR! " << __PRETTY_FUNCTION__ << " at LINE " << __LINE__ << std::endl;
	
}

template<class TBASEMAT>
void TPZMatTransientSingleSpace< TBASEMAT >::ContributeBCInterface(
    const TPZMaterialDataT<STATE> &data,
    const TPZMaterialDataT<STATE> &dataleft,
    REAL weight, 
    TPZFMatrix<STATE> &ek,
    TPZFMatrix<STATE> &ef,
    TPZBndCondT<STATE> &bc)
{
	// Mostly for implicit
	if (this->fStep == ECurrent){
		TBASEMAT::ContributeBCInterface(data,dataleft, weight,ek, ef, bc);
		return;
	}
	
	if (this->fStep == ELast){
		return;
	}
	
	
	// Mostly for explicit
	if (this->fStep == EMassMatrix){
		return;
	}
	if (this->fStep == EFluxOnly){ //Calculates ef = F-ku
		TBASEMAT::ContributeBCInterface(data, dataleft, weight,  ek, ef, bc);
		return;
	}
	
	
	PZError << "ERROR! " << __PRETTY_FUNCTION__ << " at LINE " << __LINE__ << std::endl;
	
}

template<class TBASEMAT>
void TPZMatTransientSingleSpace< TBASEMAT >::ContributeSolutionRhs(const TPZVec<STATE> &sol,
                                                             const TPZFMatrix<REAL> &phi,
                                                             REAL weight,
                                                             TPZFMatrix<STATE> &ef) {
    //Last sol is added to residual and current sol is subtracted from residual
    const STATE mult = this->fStep == ECurrent ? -1. : 1.; 
	
	const auto phr = phi.Rows();
	const auto nstate = this->NStateVariables();
	const auto DeltaT = this->TimeStep();
	for(auto i = 0; i < phr; i++) {
		for(auto k = 0; k < nstate; k++){
			ef(i*nstate+k, 0) +=
                (STATE)(mult * weight) * sol[k] * (STATE)(phi.GetVal(i,0) / DeltaT);
		}//k
	}//i
}//method

template<class TBASEMAT>
void TPZMatTransientSingleSpace< TBASEMAT >::ContributeTangent(const TPZFMatrix<REAL> &phi,
                                                         REAL weight,
                                                         TPZFMatrix<STATE> &ek) {
	const auto phr = phi.Rows();
	const auto nstate = this->NStateVariables();
	const auto DeltaT = this->TimeStep();
	for(auto i = 0; i < phr; i++) {
		for(auto j = 0; j < phr; j++){
			for(auto k = 0; k < nstate; k++){
				ek(i*nstate+k, j*nstate+k) +=
                    weight * phi.GetVal(i,0) * phi.GetVal(j,0) / DeltaT;
			}//k
		}//j
	}//i
}//method

#endif
