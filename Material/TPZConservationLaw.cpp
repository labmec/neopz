#include "TPZConsLawTest.h"
#include "pzelmat.h"
#include "pzbndcond.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzerror.h"
#include <math.h>


TPZConservationLaw::TPZConservationLaw(int nummat,REAL delta_t,int dim) : 
                                                      TPZMaterial(nummat), fDim(dim) {

  if(delta_t < 0 || delta_t > 1){
    PZError << "TPZConsLawTest::TPZConsLawTest time step parameter, default 1.0\n";
    fTimeStep = 1.0;
  }
  if(dim < 1 || dim > 3){
    PZError << "TPZConservationLaw::TPZConservationLaw (abort) error dimension = " << dim << endl;
    exit(-1);
  }
  fDim = dim;
  fTimeStep = delta_t;
}


TPZConservationLaw::~TPZConservationLaw() {
}

int TPZConservationLaw::NStateVariables() {
  return 1;
}

void TPZConservationLaw::Print(ostream &out) {
  out << "name of material : " << Name() << "\n";
  out << "properties : \n";           
  TPZMaterial::Print(out);
}

void TPZConservationLaw::Contribute(TPZVec<REAL> &x,TPZFMatrix &jacinv,TPZVec<REAL> &sol,TPZFMatrix &dsol,REAL weight,
				TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &dphi,TPZFMatrix &ek,TPZFMatrix &ef) {
  
  cout << "TPZConservationLaw::Contribute this metod does not have to be called\n";
  
}

void TPZConservationLaw::ContributeBC(TPZVec<REAL> &/*x*/,TPZVec<REAL> &/*sol*/,REAL weight,
				     TPZFMatrix &/*axes*/,TPZFMatrix &phi,TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc) {
  
  cout << "TPZConservationLaw::ContributeBC this metod does not have to be called\n";
}

/** returns the variable index associated with the name*/
int TPZConservationLaw::VariableIndex(char *name){

  cout << "TPZConservationLaw::VariableIndex this metod does not have to be called\n";
  return -1;
}

int TPZConservationLaw::NSolutionVariables(int var){
  
  if(var == 0 || var == 1 || var == 2) return 1;
  cout << "TPZConservationLaw::NSolutionVariables Error\n";
  return 0;
}

void TPZConservationLaw::Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,TPZFMatrix &/*axes*/,int var,TPZVec<REAL> &Solout){
  

}

void TPZConservationLaw::Flux(TPZVec<REAL> &/*x*/, TPZVec<REAL> &/*Sol*/, TPZFMatrix &/*DSol*/, TPZFMatrix &/*axes*/, TPZVec<REAL> &/*flux*/) {
  //Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol, TPZFMatrix &DSol, TPZFMatrix &axes, TPZVec<REAL> &flux)
}

void TPZConservationLaw::Errors(TPZVec<REAL> &/*x*/,TPZVec<REAL> &u,
			       TPZFMatrix &dudx, TPZFMatrix &axes, TPZVec<REAL> &/*flux*/,
			       TPZVec<REAL> &u_exact,TPZFMatrix &du_exact,TPZVec<REAL> &values) {
  

}
