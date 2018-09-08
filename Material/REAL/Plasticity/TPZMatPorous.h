/**
 * @file
 */

#ifndef PZPORO_H
#define PZPORO_H


#include "TPZMaterial.h"
#include "TPZMatTemporal.h" 
#include "pzporoelastoplasticmem.h"
#include "TPZMatElastoPlastic.h"
#include "pzbndcond.h"

#define TBASEPOROUS(T, TMEM) TPZMatElastoPlastic< T, TMEM >

/**
 * @brief Implements an porous media material to be used together with elastic 
 * elastoplastic mechanical counterparts.
 */
template <class T, class TMEM = TPZPoroElastoPlasticMem >
class TPZMatPorous : public TPZMatTemporal, public TPZMatElastoPlastic< T, TMEM >
{
   public:

      enum SOLUTIONVARS{ENone = -1,
		  				EPorePressure = 90};

      /** Default constructor */
      TPZMatPorous();
		
      /** Creates a material object and inserts it in the vector of
       *  material pointers of the mesh. Upon return vectorindex
       *  contains the index of the material object within the
       *  vector
       */
      TPZMatPorous(int id);

      /** Creates a material object based on the referred object and
       *  inserts it in the vector of material pointers of the mesh.
       *  Upon return vectorindex contains the index of the material
       *  object within the vector
       */
      TPZMatPorous(const TPZMatPorous<T, TMEM > &mat);
	/** @brief Default destructor */
      virtual ~TPZMatPorous();
	
      /** returns the name of the material*/
      virtual std::string Name();
	
	  /**
	   * @brief Initializes the poroelastic material coefficients
	   * @param[in] k permeability of porous medium
	   * @param[in] Mu Fluid viscosity
	   * @param[in] StorageEps poroelastic Storage Coeff at constant strain
	   * @param[in] Alpha Biot-Willis ceofficient
	   * @param[in] Rhof Fluid density
	   */
	  void SetUp(const REAL &k, const REAL &Mu, 
								const REAL &StorageEps,
								const REAL &Alpha,
								const REAL &Rhof);

      /**returns the integrable dimension of the material*/
      virtual int Dimension() const { return TBASEPOROUS(T, TMEM)::Dimension(); }

      /** returns the number of state variables associated with the material*/
      virtual int NStateVariables() { return TBASEPOROUS(T, TMEM)::NStateVariables() + 1; }

      /** print out the data associated with the material*/
      virtual void Print(std::ostream &out = std::cout, const int memory = 0);

      /**returns the variable index associated with the name*/
      virtual int VariableIndex(const std::string &name);

      /** returns the number of variables associated with the variable
	  indexed by var.  var is obtained by calling VariableIndex*/
      virtual int NSolutionVariables(int var);

      /**returns the solution associated with the var index based on
       * the finite element approximation*/
      virtual void Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout);
	
      /** Return the number of components which form the flux function
       * Method not implemented.
       */
      virtual int NFluxes() 
	  {
         PZError << "TPZMatPorous<TBASEPOROUS(T, TMEM)>::NFluxes() - Method not implemented\n";
         return 0;
      }

      /** Compute the value of the flux function to be used by ZZ error estimator.
       * Method not implemented.
       */
      virtual void Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol, TPZFMatrix<REAL> &DSol, TPZFMatrix<REAL> &axes, TPZVec<REAL> &flux)
	  {
         PZError << "TPZMatPorous<TBASEPOROUS(T, TMEM)>::Flux - Method not implemented\n";
      }

      /** Evaluate error between approximate (FEM) and exact solutions.
	   *  Method not implemented
       */
      virtual void Errors(TPZVec<REAL> &x,TPZVec<REAL> &u, TPZFMatrix<REAL> &dudx,
                          TPZFMatrix<REAL> &axes, TPZVec<REAL> &flux,
                          TPZVec<REAL> &u_exact,TPZFMatrix<REAL> &du_exact,TPZVec<REAL> &values);
      /**
       * Returns the number of norm errors: 3 (Semi H1, L2 and H1)
	   * Method not implemented
       */
      virtual int NEvalErrors() {return NStateVariables();}

      /**
       * It computes a contribution to the stiffness matrix and load vector at one integration point.
       * @param data [in] stores all input data
       * @param weight [in] is the weight of the integration rule
       * @param ek [out] is the stiffness matrix
       * @param ef [out] is the load vector
       */
      virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef);

      /**
       * It computes a contribution to the stiffness matrix and load vector at one BC integration point.
       * @param data [in] stores all input data
       * @param weight [in] is the weight of the integration rule
       * @param ek [out] is the stiffness matrix
       * @param ef [out] is the load vector
       * @param bc [in] is the boundary condition material
       */
      virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef, TPZBndCond &bc);

      /**
       * It computes a contribution to the residual vector at one integration point.
       * @param data [in] stores all input data
       * @param weight [in] is the weight of the integration rule
       * @param ef [out] is the residual vector
       */
      virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ef);

      /**
       * It computes a contribution to the stiffness matrix and load vector at one BC integration point.
       * @param data [in] stores all input data
       * @param weight [in] is the weight of the integration rule
       * @param ef [out] is the load vector
       * @param bc [in] is the boundary condition material
       */
      virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ef, TPZBndCond &bc);

      /**To create another material of the same type*/
      virtual TPZMaterial * NewMaterial();

      /** Unique identifier for serialization purposes */
      public:
virtual int ClassId() const;


      /** Save the element data to a stream */
      virtual void Write(TPZStream &buf, int withclassid) const;

      /** Read the element data from a stream */
      virtual void Read(TPZStream &buf, void *context);
	
	  /**
	   * Defining what parameters the material needs. In particular this material needs the
	   * evaluation of normal vector for the sake of boundary conditions
	   */
	  virtual void FillDataRequirements(TPZMaterialData &data);
	
	  /** Returns the porepressure and its divergent at the current integration point */
	  void ComputePorePressure(TPZMaterialData & data, REAL & Pp, TPZVec<REAL> & dPp);
	
	  /** Updates the porepressure and its divergent history */
	  void UpdatePorePressure(TPZMaterialData & data);
	
	  /** Initializes the default memory pore pressure. */
	  void SetPorePressure(const REAL Pp);
	
protected:
	
	  /** porous medium intrinsic permeability */
	  REAL fk;
	
	  /** Fluid viscosity */
	  REAL fMu;
	
	  /** Storage coefficient at constant strain */
	  REAL fStorageEps;
	
	  /** Biot-Willis poroelastic coefficient */
	  REAL fAlpha;
	
	  /** fluid density */
	  REAL fRhof;
	
};

template <class T, class TMEM>
int TPZMatPorous<T, TMEM >::ClassId() const{
    return Hash("TPZMatPorous") ^ TPZMatTemporal::ClassId() << 1 ^ TPZMatElastoPlastic<T, TMEM>::ClassId() << 2;
}
#endif
