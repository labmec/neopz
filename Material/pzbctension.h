// -*- c++ -*-
// $Id: pzbctension.h,v 1.3 2003-11-04 00:26:33 phil Exp $

#ifndef BCTENSIONHPP
#define BCTENSIONHPP

#include "pzbndcond.h"
#include "TPZPlacaOrthotropic.h"

using namespace std;

template <class T, int N>
class TPZManVector;
class TPZInterpolatedElement;
class TPZMulticamadaOrthotropic;


// this class defines the boundary condition for a 1d linear problem

class TPZBCTension : public TPZBndCond {

  REAL fHight,fLxLyDim,fMx,fMy,fMxy,fNx,fNy,fNxy,fQx,fQy,fQxy;
  TPZMulticamadaOrthotropic *fMultCam;

  private:

  TPZInterpolatedElement *fReference;

  public :
    
    ~TPZBCTension(){}

  TPZBCTension(TPZMaterial *material,int id,int type,TPZFMatrix &val1,TPZFMatrix &val2, TPZMulticamadaOrthotropic *mult);


  void FeedBCTension(){;}

  void Tensor(TPZVec<REAL> &Point, int placa, int var, TPZManVector<REAL> &SolOut);


  virtual int NFluxes(){ return Material()->NFluxes(); }

  int NStateVariables() { return Material()->NStateVariables(); }

  void Contribute(TPZVec<REAL> &x,TPZFMatrix &,TPZVec<REAL> &sol,TPZFMatrix &dsol, REAL
		  weight,TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &dphi,TPZFMatrix &ek,TPZFMatrix &ef) {
    if(fForcingFunction) {
      TPZManVector<REAL> result(Val2().Rows());
      fForcingFunction(x,result);
      int i;
      for(i=0; i<Val2().Rows(); i++) {
	Val2()(i,0) = result[i];
      }
    }


    Material()->ContributeBC(x,sol,weight,axes,phi,ek,ef,*this);

  }

  void ContributeBC(TPZVec<REAL> &x,TPZVec<REAL> &sol,REAL weight,TPZFMatrix &axes,
		    TPZFMatrix &phi,TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc) {
  }
};


#endif


