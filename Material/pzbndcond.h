// -*- c++ -*-

//$Id: pzbndcond.h,v 1.8 2004-04-06 14:55:43 erick Exp $

//HEADER FILE FOR CLASS BNDCOND

#ifndef BNDCONDHPP
#define BNDCONDHPP


#include <iostream>

using namespace std;

#include "pzreal.h"
#include "pzdiscgal.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include "pzfmatrix.h"

#ifdef _AUTODIFF
#include "fadType.h"
#endif

template <class T, int N>
class TPZManVector;

// this class defines the boundary condition for a 1d linear problem

class TPZBndCond : public TPZDiscontinuousGalerkin {

  friend class TPZMaterial;
protected:
  int 		fType;		              // boundary condition type
  TPZFMatrix	fBCVal1;            // first value of boundary condition
  TPZFMatrix	fBCVal2;            // second value of boundary condition
  TPZMaterial	*fMaterial;	        // pointer to material which created bc

  TPZBndCond(TPZBndCond & bc) : TPZDiscontinuousGalerkin(bc), fBCVal1(bc.fBCVal1),
    fBCVal2(bc.fBCVal2){
    fMaterial = bc.fMaterial;
    fType = bc.fType;
  }

  public :

    ~TPZBndCond(){}

  TPZBndCond(TPZMaterial *material,int id,int type,TPZFMatrix &val1,TPZFMatrix &val2) :
    TPZDiscontinuousGalerkin(id), fBCVal1(val1), fBCVal2(val2) {
    //cria um novo material
    fMaterial = material;
    fType = type;

  }

  TPZBndCond(TPZBndCond &copy, TPZMaterial *ref) : TPZDiscontinuousGalerkin(copy), fType(copy.fType),
						   fBCVal1(copy.fBCVal1), fBCVal2(copy.fBCVal2), fMaterial(ref) {}
 
  void SetMaterial(TPZMaterial * mat) { fMaterial = mat;}

  /**returns the integrable dimension of the material*/
  int Dimension() { return fMaterial->Dimension(); }

  virtual int NFluxes(){ return fMaterial->NFluxes(); }

  int NStateVariables() { return fMaterial->NStateVariables(); }

  /**
   * Returns the number of norm errors. Default is 3: energy, L2 and H1.
   */
  virtual int NEvalErrors() {return fMaterial->NEvalErrors();}

  int Type() { return fType; }

  TPZFMatrix &Val1() { return fBCVal1; }

  TPZFMatrix &Val2() { return fBCVal2; }

  TPZMaterial *Material() { return fMaterial; }

  /**compute the value of the flux function to be used by ZZ error estimator*/
  void Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol, TPZFMatrix &DSol, TPZFMatrix &axes, TPZVec<REAL> &flux){
    flux.Fill(0.);
  }

  void Print(ostream & out = cout) {
    out << " Boundary condition number = " << Id() << "\n";
    out << " boundary condition type = " << fType << "\n";
    out << " val1 = \n"; fBCVal1.Print("fBCVal1",out);
    out << " val2 = \n"; fBCVal2.Print("fBCVal2",out);
  }

  void Contribute(TPZVec<REAL> &x,TPZFMatrix &,TPZVec<REAL> &sol,TPZFMatrix &dsol, REAL
	weight,TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &dphi,TPZFMatrix &ek,TPZFMatrix &ef) {
    if(fForcingFunction) {
      TPZManVector<REAL> result(fBCVal2.Rows());
      fForcingFunction(x,result);
      int i;
      for(i=0; i<fBCVal2.Rows(); i++) {
	fBCVal2(i,0) = result[i];
      }
    }
    //clone meshes required analysis
    int typetmp = fType;
    if (fType == 50){
	    int i;
	    for (i=0;i<sol.NElements();i++){
		fBCVal2(i,0) = gBigNumber*sol[i];
		fBCVal1(i,i) = gBigNumber;
	    }
	    fType = 2;
    }
    fMaterial->ContributeBC(x,sol,weight,axes,phi,ek,ef,*this);
    fType = typetmp;
  }


#ifdef _AUTODIFF

  void ContributeEnergy(TPZVec<REAL> &x,
			      TPZVec<FADFADREAL> &sol,
			      TPZVec<FADFADREAL> &dsol,
			      FADFADREAL &U,
			      REAL weight);

#endif


  void ContributeBC(TPZVec<REAL> &x,TPZVec<REAL> &sol,REAL weight,TPZFMatrix &axes,
		    TPZFMatrix &phi,TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc) {
  }


  void Errors(TPZVec<REAL> &x,TPZVec<REAL> &sol,TPZFMatrix &dsol, TPZFMatrix &axes, TPZVec<REAL> &flux,
	      TPZVec<REAL> &uexact,TPZFMatrix &duexact,TPZVec<REAL> &val){
    val.Fill(0.);
  }

  virtual void Clone(TPZAdmChunkVector<TPZMaterial *> &matvec);

  virtual void ContributeInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL,TPZVec<REAL> &solR,TPZFMatrix &dsolL,
				   TPZFMatrix &dsolR,REAL weight,TPZVec<REAL> &normal,TPZFMatrix &phiL,
				   TPZFMatrix &phiR,TPZFMatrix &dphiL,TPZFMatrix &dphiR,
				   TPZFMatrix &ek,TPZFMatrix &ef);

  virtual void ContributeInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL,TPZVec<REAL> &solR,TPZFMatrix &dsolL,
				   TPZFMatrix &dsolR,REAL weight,TPZVec<REAL> &normal,TPZFMatrix &phiL,
				   TPZFMatrix &phiR,TPZFMatrix &dphiL,TPZFMatrix &dphiR,
				   TPZFMatrix &ef);
  
  virtual void ContributeInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL,TPZVec<REAL> &solR,TPZFMatrix &dsolL,
				   TPZFMatrix &dsolR,REAL weight,TPZVec<REAL> &normal,TPZFMatrix &phiL,
				   TPZFMatrix &phiR,TPZFMatrix &dphiL,TPZFMatrix &dphiR,
				   TPZFMatrix &ek,TPZFMatrix &ef, int LeftPOrder, int RightPOrder, REAL faceSize);

  virtual void ContributeBCInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL, TPZFMatrix &dsolL, REAL weight, TPZVec<REAL> &normal,
    TPZFMatrix &phiL,TPZFMatrix &dphiL, TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc) {
    //NOTHING TO BE DONE HERE
  }

};


#endif
