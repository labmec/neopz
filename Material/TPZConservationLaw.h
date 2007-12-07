#ifndef CONSERVATIONLAWHPP
#define CONSERVATIONLAWHPP

#include <iostream>
#include "pzmaterial.h"
#include "pzfmatrix.h"
#include "pzvec.h"

/// This class implements the interface for conservation laws, keeping track of the timestep as well
/**
THIS CLASS IS DEPRECATED BY TPZConservationLaw2
*/
class TPZConservationLaw  : public TPZMaterial {

  int fDim;
  REAL fTimeStep;

  public :

  REAL fDelta;

  TPZConservationLaw(int nummat,REAL delta_t,int dim);

  /**copy constructor*/
  TPZConservationLaw(TPZConservationLaw &copy);

  /**To create another material of the same type*/
  TPZAutoPointer<TPZMaterial> NewMaterial();

  ~TPZConservationLaw(){};

  /**
   * compute the boundary condition left solution
   */
  virtual void ComputeSolLeft(TPZVec<REAL> &solr,TPZVec<REAL> &soll,TPZVec<REAL> &normal,TPZBndCond *bcleft);

  /**
   * compute the boundary condition right solution
   */
  virtual void ComputeSolRight(TPZVec<REAL> &solr,TPZVec<REAL> &soll,TPZVec<REAL> &normal,TPZBndCond *bcright);

  /**
   * termodinamic pressure determined by the law of an ideal gas
   */
  virtual REAL Pressure(TPZVec<REAL> &U);

  virtual REAL Gamma();

  //virtual int IntegrationDegree() = 0;

  //virtual void SetIntegDegree(int degree) = 0;

  virtual void SetDelta(REAL delta){fDelta = delta;}

  virtual void SetDeltaTime(REAL maxveloc,REAL deltax,int degree);

  REAL Delta();

  virtual void SetTimeStep(REAL timestep){fTimeStep = timestep;}

  virtual REAL TimeStep(){return fTimeStep;}

  int Dimension() { return fDim;}

  int NStateVariables();

  virtual void Print(std::ostream & out);

  virtual std::string Name() { return "TPZConservationLaw"; }

  virtual void Contribute(TPZMaterialData &data,
                            REAL weight,
                            TPZFMatrix &ek,
                            TPZFMatrix &ef);

  virtual void ContributeInterface(TPZMaterialData &data,
                                     REAL weight,
                                     TPZFMatrix &ek,
                                     TPZFMatrix &ef);

  virtual void ContributeBC(TPZMaterialData &data,
                              REAL weight,
                              TPZFMatrix &ek,
                              TPZFMatrix &ef,
                              TPZBndCond &bc);

  virtual void Contribute(TPZMaterialData &data,
                            REAL weight,
                            TPZFMatrix &ef);

  virtual int VariableIndex(char *name);

  virtual int NSolutionVariables(int var);

  virtual int NFluxes(){ return 1;}

  virtual void Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,TPZFMatrix &axes,int var,TPZVec<REAL> &Solout);

  /**compute the value of the flux function to be used by ZZ error estimator*/
  virtual void Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol, TPZFMatrix &DSol, TPZFMatrix &axes, TPZVec<REAL> &flux);

  void Errors(TPZVec<REAL> &x,TPZVec<REAL> &u,
	      TPZFMatrix &dudx, TPZFMatrix &axes, TPZVec<REAL> &flux,
	      TPZVec<REAL> &u_exact,TPZFMatrix &du_exact,TPZVec<REAL> &values);
};


// inline void TPZConservationLaw::ContributeInterface(TPZMaterialData &data,
//                                                       REAL ,
//                                                       TPZFMatrix &){
//   PZError << "TPZConservationLaw::ContributeInterface it would never have to be called\n";
// }

inline void TPZConservationLaw::Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,TPZFMatrix &axes,int var,TPZVec<REAL> &Solout){
  PZError << "TPZConservationLaw::Solution it would never have to be called\n";
}

inline void TPZConservationLaw::SetDeltaTime(REAL maxveloc,REAL deltax,int degree){
  PZError << "TPZConservationLaw::SetDeltaTime it would never have to be called\n";
}

inline REAL TPZConservationLaw::Pressure(TPZVec<REAL> &U) {
  PZError << "TPZConservationLaw::Pressure it would never have to be called\n";
  return 0.0;
}

inline REAL TPZConservationLaw::Gamma() {
  PZError << "TPZConservationLaw::Gamma it would never have to be called\n";
  return 0.0;
}

inline void TPZConservationLaw::ComputeSolLeft(TPZVec<REAL> &,TPZVec<REAL> &,TPZVec<REAL> &,TPZBndCond *){
  PZError << "TPZConservationLaw::ComputeSolLeft it would never have to be called\n";
}

inline void TPZConservationLaw::ComputeSolRight(TPZVec<REAL> &,TPZVec<REAL> &,TPZVec<REAL> &,TPZBndCond *){
  PZError << "TPZConservationLaw::ComputeSolRight it would never have to be called\n";
}
#endif
