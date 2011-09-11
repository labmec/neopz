/**
 * \file
 * @brief DEPRECATED FILE. Contains implementations of the TPZConservationLaw methods.
 */
#include "TPZConservationLaw.h"
#include "pzelmat.h"
#include "pzbndcond.h"
#include "pzmatrix.h" 
#include "pzfmatrix.h"
#include "pzerror.h"
#include <math.h>
#include <stdlib.h>

using namespace std;

TPZConservationLawDEP::TPZConservationLawDEP(int nummat,REAL delta_t,int dim) : 
TPZMaterial(nummat), fDim(dim) {
	
	fTimeStep = delta_t;
	if(delta_t < 0 || delta_t > 1){
		PZError << "TPZConservationLawDEP::TPZConservationLawDEP time step parameter > 1 , default 1.0\n";
		fTimeStep = 1.0;
	}
	if(dim < 1 || dim > 3){
		PZError << "TPZConservationLawDEP::TPZConservationLawDEP (abort) error dimension = " << dim << endl;
		exit(-1);
	}
	fDim = dim;
}

TPZConservationLawDEP::TPZConservationLawDEP(TPZConservationLawDEP &copy) : TPZMaterial(copy) {
	fDim = copy.fDim;
	fTimeStep = copy.fTimeStep;
	fDelta = copy.fDelta;
}

TPZAutoPointer<TPZMaterial> TPZConservationLawDEP::NewMaterial() {
	PZError << "TPZConservationLawDEP::::NewMaterial is called.\n";
	return 0;
}

int TPZConservationLawDEP::NStateVariables() {
	return 1;
}

void TPZConservationLawDEP::Print(std::ostream &out) {
	out << "name of material : " << Name() << "\n";
	out << "properties : \n";           
	TPZMaterial::Print(out);
}

void TPZConservationLawDEP::Contribute(TPZMaterialData &data,
                                    REAL weight,
                                    TPZFMatrix &ek,
                                    TPZFMatrix &ef) {
	
	cout << "TPZConservationLawDEP::Contribute this metod does not have to be called\n";
	
}

void TPZConservationLawDEP::Contribute(TPZMaterialData &data,
                                    REAL weight,
                                    TPZFMatrix &ef) {
	
	cout << "TPZConservationLawDEP::Contribute this metod does not have to be called\n";
	
}

void TPZConservationLawDEP::ContributeBC(TPZMaterialData &data,
                                      REAL weight,
                                      TPZFMatrix &ek,
                                      TPZFMatrix &ef,
                                      TPZBndCond &bc) {
	
	cout << "TPZConservationLawDEP::ContributeBC this metod does not have to be called\n";
}

void TPZConservationLawDEP::ContributeInterface(TPZMaterialData &data,
											 REAL weight, 
											 TPZFMatrix &ek,
											 TPZFMatrix &ef)
{
	cout << __PRETTY_FUNCTION__ << " this metod should not be called\n";
}
/** returns the variable index associated with the name*/
int TPZConservationLawDEP::VariableIndex(const std::string &name){
	
	cout << "TPZConservationLawDEP::VariableIndex this metod does not have to be called\n";
	return -1;
}

int TPZConservationLawDEP::NSolutionVariables(int var){
	
	if(var == 0 || var == 1 || var == 2) return 1;
	cout << "TPZConservationLawDEP::NSolutionVariables Error\n";
	return 0;
}

// void TPZConservationLawDEP::Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,TPZFMatrix &/*axes*/,int var,TPZVec<REAL> &Solout){

//   cout << "TPZConservationLawDEP::Solution nao deve ser chamada\n";

// }

void TPZConservationLawDEP::Flux(TPZVec<REAL> &/*x*/, TPZVec<REAL> &/*Sol*/, TPZFMatrix &/*DSol*/, TPZFMatrix &/*axes*/, TPZVec<REAL> &/*flux*/) {
	//Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol, TPZFMatrix &DSol, TPZFMatrix &axes, TPZVec<REAL> &flux)
}

void TPZConservationLawDEP::Errors(TPZVec<REAL> &/*x*/,TPZVec<REAL> &u,
								TPZFMatrix &dudx, TPZFMatrix &axes, TPZVec<REAL> &/*flux*/,
								TPZVec<REAL> &u_exact,TPZFMatrix &du_exact,TPZVec<REAL> &values) {
	
	
}
