/** @file pzmaterial.h
 *
 * Header file for abstract class TPZMaterial.
 */

#ifndef PZMATERIALHPP
#define PZMATERIALHPP

#include <iostream>
#include <string.h>

#include "pzreal.h"
#include "pzvec.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzadmchunk.h"

#ifdef _AUTODIFF
#include "fadType.h"
#endif

using namespace std;

class TPZBndCond;
class TPZMaterial;

class  TPZMaterial
{
      int fId;

   protected:
      void (*fForcingFunction)(TPZVec<REAL> &loc,TPZVec<REAL> &result);

   public:

      static REAL gBigNumber;

      /** Creates a material object and inserts it in the vector of
       *  material pointers of the mesh. Upon return vectorindex
       *  contains the index of the material object within the
       *  vector
       */
      TPZMaterial(int id);

      /** Creates a material object based on the referred object and
       *  inserts it in the vector of material pointers of the mesh.
       *  Upon return vectorindex contains the index of the material
       *  object within the vector
       */
      TPZMaterial(TPZMaterial &mat);

      virtual ~TPZMaterial() {}

      /** returns the name of the material*/
      virtual char *Name() { return "no_name"; }

      /**returns the integrable dimension of the material*/
      virtual int Dimension() = 0;

      int Id() const { return fId; }
      void SetId(int id) { fId = id; }

      /** returns the number of state variables associated with the material*/
      virtual int NStateVariables() = 0;

      /** return the number of components which form the flux function*/
      virtual int NFluxes() {return 0;}

      /** print out the data associated with the material*/
      virtual void Print(ostream &out = std::cout);

      /**returns the variable index associated with the name*/
      virtual int VariableIndex(char *name);

      /** returns the number of variables associated with the variable
	  indexed by var.  var is obtained by calling VariableIndex*/
      virtual int NSolutionVariables(int var);

      /**returns the solution associated with the var index based on
       * the finite element approximation*/
      virtual void Solution(TPZVec<REAL> &Sol, TPZFMatrix &DSol,
			    TPZFMatrix &axes, int var, TPZVec<REAL> &Solout);

      /**compute the value of the flux function to be used by ZZ error
       * estimator*/
      virtual void Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol, 
			TPZFMatrix &DSol, TPZFMatrix &axes,
			TPZVec<REAL> &flux) {}

      /**Create an object TPZBndCond derived of TPZMaterial*/
      virtual TPZBndCond *CreateBC(int id, int typ, TPZFMatrix &val1,
				   TPZFMatrix &val2);

      /**Compute contribution to the stiffness matrix and right hand
       * side at an integration point*/
      virtual void Contribute(TPZVec<REAL> &x, TPZFMatrix &jacinv,
			      TPZVec<REAL> &sol, TPZFMatrix &dsol,
			      REAL weight,TPZFMatrix &axes,
			      TPZFMatrix &phi, TPZFMatrix &dphi,
			      TPZFMatrix &ek, TPZFMatrix &ef) = 0;

#ifdef _AUTODIFF

      /**Compute contribution to the energy at an integration point*/
      virtual void ContributeEnergy(TPZVec<REAL> &x,
			      TPZVec<FADFADREAL> &sol,
			      TPZVec<FADFADREAL> &dsol,
			      FADFADREAL &U,
			      REAL weight);

#endif


      /** Compute contribution to the stiffness matrix and right hand
       * side at the integration point of a boundary*/
      virtual void ContributeBC(TPZVec<REAL> &x, TPZVec<REAL> &sol,
				REAL weight, TPZFMatrix &axes,
				TPZFMatrix &phi, TPZFMatrix &ek,
				TPZFMatrix &ef, TPZBndCond &bc) = 0;

#ifdef _AUTODIFF

      /** Compute contribution of BC to the Energy*/
      virtual void ContributeBCEnergy(TPZVec<REAL> & x,
	TPZVec<FADFADREAL> & sol, FADFADREAL &U,
	REAL weight, TPZBndCond &bc);

#endif
      /**Set a procedure as source function for the material.  loc
	 corresponds to the coordinate of the point where the source
	 function is applied*/
      void SetForcingFunction(void (*fp)(TPZVec<REAL> &loc,
					 TPZVec<REAL> &result))
      {
	 fForcingFunction = fp;
      }

      /**Compute the error due to the difference between the
	 interpolated flux and the flux computed based on the
	 derivative of the solution*/
      virtual void Errors(TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix &dsol,
			  TPZFMatrix &axes, TPZVec<REAL> &flux,
			  TPZVec<REAL> &uexact, TPZFMatrix &duexact,
			  TPZVec<REAL> &val) {}

      /**To create another material of the same type*/
      virtual TPZMaterial *NewMaterial();

      /**Read data of the material from a istream (file data)*/
      virtual void SetData(istream &data);

      /**Compute contribution to the stiffness matrix and right hand
       * side at an integration point*/
      virtual void Contribute(TPZVec<REAL> &x, TPZFMatrix &jacinv,
			      TPZVec<REAL> &sol, TPZFMatrix &dsol, REAL weight,
			      TPZFMatrix &axes, TPZFMatrix &phi,
			      TPZFMatrix &dphi, TPZFMatrix &ef);

      /**
       * Create a copy of the material object and put it in the vector
       * which is passed on
       */
      virtual void Clone(TPZAdmChunkVector<TPZMaterial *> &matvec);

      /**To return a numerical flux type to apply over the interfaces
       * of the elements*/
      virtual int FluxType() { return 2; }

      /** Factor to diffussive term*/
      virtual int IdBC(double *x) { return 5; }

      /**Compute the contribution to stiffness matrix and load vector
	 of the fluxes over the interface between two computational
	 elements*/
      virtual void ContributeOverInterface(TPZVec<REAL> &x,
					   TPZVec<REAL> &leftsol,
					   TPZVec<REAL> &rightsol,
					   REAL weight, REAL area, int type,
					   TPZFMatrix &axes,
					   TPZFMatrix &phileft,
					   TPZFMatrix &phiright,
					   TPZVec<REAL> &normal,
					   TPZFMatrix &ek, TPZFMatrix &ef);
};

inline void TPZMaterial::ContributeOverInterface(
   TPZVec<REAL> &/*x*/, TPZVec<REAL> &/*leftsol*/, TPZVec<REAL> &/*rightsol*/,
   REAL /*weight*/, REAL /*area*/, int /*type*/, TPZFMatrix &/*axes*/,
   TPZFMatrix &/*phileft*/, TPZFMatrix &/*phiright*/, TPZVec<REAL> &/*normal*/,
   TPZFMatrix &/*ek*/, TPZFMatrix &/*ef*/)
{
   cout << "TPZMaterial::ContributeOverInterface is called.\n";
}

#endif

