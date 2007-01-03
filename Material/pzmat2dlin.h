#ifndef MAT2DLINHPP
#define MAT2DLINHPP

#include <iostream>
#include "pzmaterial.h"
#include "pzfmatrix.h"
//#include "pzreal.h"
struct TPZElementMatrix;
class TPZBndCond;
template<class T>
class TPZVec;

//const Float BIGNUMBER = 1.e9;

class TPZMat2dLin : public TPZMaterial{

  TPZFMatrix    fKxx, fKxy, fKyx, fKyy, fKx0, fK0x, fKy0, fK0y, fK00, fXf;
  public :

    TPZMat2dLin(int num = 0) : TPZMaterial(num), fKxx(), fKxy(),
     fKyx() , fKyy(), fKx0(), fK0x(), fKy0(), fK0y(), fK00(), fXf() {
    }

  TPZMat2dLin(TPZMat2dLin &copy) : TPZMaterial(copy), 
	  fKxx(copy.fKxx), fKxy(copy.fKxy), fKyx(copy.fKyx), fKyy(copy.fKyy), 
	  fKx0(copy.fKx0), fK0x(copy.fK0x), fKy0(copy.fKy0), 
	  fK0y(copy.fK0y), fK00(copy.fK00), fXf(copy.fXf){}

  virtual int NStateVariables() { return fKxx.Rows(); }

//  int NFluxes() { return NStateVariables(); }

  int Dimension() { return 2; }

  void Print(std::ostream & out = std::cout);

  void SetMaterial(TPZFMatrix &xkin,TPZFMatrix &xcin,TPZFMatrix &xfin){
	  int r = xkin.Rows();
    fKxx = xkin;
	fKyy = xkin;
	fK00 = xcin;
    fXf = xfin;
	fKxy.Redim(r,r);
	fKyx.Redim(r,r);
	fKx0.Redim(r,r);
	fK0x.Redim(r,r);
	fKy0.Redim(r,r);
	fK0y.Redim(r,r);
  }

  void ConvectionDiffusion(REAL angle,REAL diff);

  TPZFMatrix &Xk() {return fKxx;}
  TPZFMatrix &Ck() {return fK00;}
  TPZFMatrix &Xf() {return fXf;}  

  virtual char *Name() { return "TPZMat2dLin"; }

//  virtual TPZBndCond *CreateBC(int num,int typ,TPZFMatrix &val1,TPZFMatrix &val2);

  virtual void Contribute(TPZVec<REAL> &x,TPZFMatrix &,TPZVec<REAL> &sol,TPZFMatrix &,REAL weight,
			  TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &dphi,TPZFMatrix &ek,TPZFMatrix &ef);

  virtual void ContributeBC(TPZVec<REAL> &x,TPZVec<REAL> &sol,REAL weight,
			    TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc);

  virtual int NFluxes();

  virtual void Flux(TPZVec<REAL> &x, TPZVec<REAL> &u, TPZFMatrix &dudx, TPZFMatrix &axes, TPZVec<REAL> &fl);

  void Errors(TPZVec<REAL> &x,TPZVec<REAL> &u,TPZFMatrix &dudx,TPZFMatrix &axes,TPZVec<REAL> &flux,
			 TPZVec<REAL> &u_exact,TPZFMatrix &du_exact,TPZVec<REAL> &values);

  virtual int VariableIndex(char *name);

  virtual int NSolutionVariables(int index);

  void Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,TPZFMatrix &axes, int var,TPZVec<REAL> &Solout);

  /**
   * Create a copy of the material object
   */
  virtual TPZMaterial *NewMaterial();

TPZBndCond *OutflowFlux(TPZAutoPointer<TPZMaterial> &reference, int bc);

  /**
  * returns the unique identifier for reading/writing objects to streams
  */
  virtual int ClassId() const;
  /**
  Save the element data to a stream
  */
  virtual void Write(TPZStream &buf, int withclassid);
  
  /**
  Read the element data from a stream
  */
  virtual void Read(TPZStream &buf, void *context);


};

#endif


