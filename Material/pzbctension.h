// -*- c++ -*-
// $Id: pzbctension.h,v 1.6 2003-11-07 00:40:34 phil Exp $

#ifndef BCTENSIONHPP
#define BCTENSIONHPP

#include "pzbndcond.h"
#include "TPZPlacaOrthotropic.h"
#include "TPZMulticamadaOrtho.h"

using namespace std;

template <class T, int N>
class TPZManVector;
class TPZInterpolatedElement;
class TPZMulticamadaOrthotropic;


// this class defines the boundary condition for a 1d linear problem

class TPZBCTension : public TPZBndCond {

  TPZMulticamadaOrthotropic *fMultCam;
  int fCamada;
  REAL fSign;

  private:

  //  TPZInterpolatedElement *fReference;

  public :
    
    ~TPZBCTension(){}

  TPZBCTension(TPZMaterial *material,int id,int type,TPZFMatrix &val1,TPZFMatrix &val2, REAL sign, TPZMulticamadaOrthotropic *mult, int camada);



  virtual int NFluxes(){ return Material()->NFluxes(); }

  int NStateVariables() { return Material()->NStateVariables(); }

  void Contribute(TPZVec<REAL> &x,TPZFMatrix &jacinv,TPZVec<REAL> &sol,TPZFMatrix &dsol, REAL
		  weight,TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &dphi,TPZFMatrix &ek,TPZFMatrix &ef) {

    int typekeep = fType;
    if(fType == 4) {
      TPZManVector<REAL,3> normal(3);
      normal[0] = fSign*(axes(0,1)*axes(1,2)-axes(0,2)*axes(1,1));
      normal[1] = fSign*(axes(0,2)*axes(1,0)-axes(0,0)*axes(1,2));
      normal[2] = fSign*(axes(0,0)*axes(1,1)-axes(0,1)*axes(1,0));
      TPZFNMatrix<9> tensor(3,3);
      fMultCam->Tensor(x,fCamada,tensor);
      int i,j;
      Val2().Zero();
      for(i=0; i<3; i++) {
	for(j=0; j<3; j++) {
	  Val2()(i,0) += tensor(i,j)*normal[j];
	}
      }
      fType = 1;
      Material()->ContributeBC(x,sol,weight,axes,phi,ek,ef,*this);
      fType = typekeep;
    } else {
      TPZBndCond::Contribute(x,jacinv,sol,dsol,weight,axes,phi,dphi,ek,ef);
    }

  }
};


#endif


