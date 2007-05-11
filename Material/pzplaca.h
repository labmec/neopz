#ifndef PLACAHPP
#define PLACAHPP

#include "pzmaterial.h"
#include "pzfmatrix.h"
#include "pzvec.h"
//#include "pzreal.h"
class TPZBndCond;
template<class T>
class TPZVec;

//const Float BIGNUMBER = 1.e9;

class TPZPlaca : public TPZMaterial{

  TPZFMatrix fnaxes;
  REAL fE1,fE2,fG12,fG13,fG23,fh,ff,fmi,fni1,fni2;
  TPZFMatrix fRmat, fRmatT;
  TPZFMatrix fKxxR,fKyyR, fKxyR, fKyxR, fBx0R, fB0xR,
             fBy0R,fB0yR, fB00R;
  TPZVec<REAL> fXF;
  public :

  TPZPlaca(int num, REAL h, REAL f, REAL E1 , REAL E2 ,
           REAL ni1 , REAL ni2 , REAL G12 , REAL G13 ,
           REAL G23 , TPZFMatrix &naxes, TPZVec<REAL> &xf);

  virtual int NStateVariables() { return 6; }

//  int NFluxes() { return NStateVariables(); }

  int Dimension() { return 2; }

  void Print(std::ostream & out);

  virtual char *Name() { return "TPZPlaca"; }

//  virtual TPZBndCond *CreateBC(int num,int typ,TPZFMatrix &val1,TPZFMatrix &val2);

  virtual void Contribute(TPZMaterialData &data,
                            REAL weight,
                            TPZFMatrix &ek,
                            TPZFMatrix &ef);

  virtual void ContributeBC(TPZMaterialData &data,
                              REAL weight,
                              TPZFMatrix &ek,
                              TPZFMatrix &ef,
                              TPZBndCond &bc);

  virtual int NFluxes();

  virtual void Flux(TPZVec<REAL> &x, TPZVec<REAL> &u, TPZFMatrix &dudx, TPZFMatrix &axes, TPZVec<REAL> &fl);

  virtual void Errors(TPZVec<REAL> &x,TPZVec<REAL> &u,TPZFMatrix &dudx, TPZFMatrix &axes, TPZVec<REAL> &flux,
		      TPZVec<REAL> &u_exact,TPZFMatrix &du_exact,TPZVec<REAL> &values);

  /**returns the variable index associated with the name*/
  virtual int VariableIndex(char *name);

  /** returns the number of variables associated with the variable indexed by var.
      var is obtained by calling VariableIndex*/
  virtual int NSolutionVariables(int var);

  /**returns the solution associated with the var index based on the finite element approximation*/
  virtual void Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,TPZFMatrix &axes,int var,TPZVec<REAL> &Solout);

  /**Exact solution for tests*/
  void SetExactFunction( void (*fp)(TPZFMatrix &axes,TPZVec<REAL> &x,TPZFMatrix &uexact,TPZFMatrix &duexact) )
    {
      fExactFunction = fp;
    }
  
 protected:

  void (*fExactFunction)(TPZFMatrix &axes,TPZVec<REAL> &x,TPZFMatrix &uexact,TPZFMatrix &duexact);

};

#endif
