/**
 * @file
 */

#ifndef PZPOSTPROCMAT_H
#define PZPOSTPROCMAT_H

#include "TPZMaterial.h"

#include "tpzautopointer.h"
#include "pzvec.h"

#include <string>

  /**
   * @brief Implements an elastoplastic material and uses the memory feature to store the damage variables
   * This material works only together with the Plasticity Library.
   */
class TPZPostProcVar : public TPZSavable
{
public:
		
	  TPZPostProcVar():fIndex(-1), fName(""), fNumEq(-1) {}

	  TPZPostProcVar(const TPZPostProcVar & source)
	  {
	      this->operator=(source);
	  }
	
	  TPZPostProcVar& operator=(const TPZPostProcVar & source)
	  {
		  fIndex = source.fIndex;
		  fName = source.fName;
		  fNumEq  = source.fNumEq;
		  return *this;
	  }
		
	  ~TPZPostProcVar(){}	int ClassId() const override;
        void Read(TPZStream &buf, void *context) override;
        void Write(TPZStream &buf, int withclassid) const override;

public: //members
	
	  int64_t fIndex;
	
	  std::string fName;
	
	  int64_t fNumEq;
};

class  TPZPostProcMat : public TPZMaterial
{
   public:
		
      /** @brief Default constructor */
      TPZPostProcMat();		
		
      /**
	   * @brief Creates a material object and inserts it in the vector of
       *  material pointers of the mesh. Upon return vectorindex
       *  contains the index of the material object within the
       *  vector
       */
      TPZPostProcMat(int64_t id);

      /**
	   * @brief Creates a material object based on the referred object and
       *  inserts it in the vector of material pointers of the mesh.
       */
      TPZPostProcMat(const TPZPostProcMat &mat);

      virtual ~TPZPostProcMat();
	
      /** @brief returns the name of the material*/
      virtual std::string Name() override;

      /** @brief returns the integrable dimension of the material*/
      virtual int Dimension() const override;

      /** @brief returns the number of state variables associated with the material*/
      virtual int NStateVariables() const override;

      /** @brief print out the data associated with the material*/
      virtual void Print(std::ostream &out = std::cout) override;

      /** @brief returns the variable index associated with the name*/
      virtual int VariableIndex(const std::string &name) override;

      /**
	   * @brief returns the number of variables associated with the variable
	   * indexed by var.  var is obtained by calling VariableIndex
	   */
      virtual int NSolutionVariables(int var) override;

      /**
	   * @brief returns the solution associated with the var index based on
       * the finite element approximation
	   */
      virtual void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout) override;
	
      /**
       * @brief It computes a contribution to the stiffness matrix and load vector at one integration point.
       * @param data [in] stores all input data
       * @param weight [in] is the weight of the integration rule
       * @param ek [out] is the stiffness matrix
       * @param ef [out] is the load vector
       */
      virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;

      /**
       * @brief It computes a contribution to the residual vector at one integration point.
       * @param data [in] stores all input data
       * @param weight [in] is the weight of the integration rule
       * @param ef [out] is the residual vector
       */
      virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef) override;


      /**
       * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
       * @param data [in] stores all input data
       * @param weight [in] is the weight of the integration rule
       * @param ek [out] is the stiffness matrix
       * @param ef [out] is the load vector
       * @param bc [in] is the boundary condition material
       * @since April 16, 2007
       */
      virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) override;
	
      /**
       * @brief It computes a contribution to stiffness matrix and load vector at one integration point
       * @param data [in]
       * @param dataleft [in]
       * @param dataright [in]
       * @param weight [in]
       * @param ek [out] is the stiffness matrix
       * @param ef [out] is the load vector
       */
      virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, 
                                       REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;

      /**
       * @brief It computes a contribution to stiffness matrix and load vector at one BC integration point
       * @param data [in]
       * @param dataleft [in]
       * @param weight [in]
       * @param ek [out] is the stiffness matrix
       * @param ef [out] is the load vector
       * @param bc [in] is the boundary condition object
       */
      virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, 
                                         REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override;

      /** @brief Unique identifier for serialization purposes */
      public:
        virtual int ClassId() const override;


      /** @brief Save the element data to a stream */
      virtual void Write(TPZStream &buf, int withclassid) const override;

      /** @brief Read the element data from a stream */
      virtual void Read(TPZStream &buf, void *context) override;
	
	  /**
	   * @brief Defining what parameters the material needs. In particular this material needs the
	   * evaluation of normal vector for the sake of boundary conditions
	   */
	  virtual void FillDataRequirements(TPZMaterialData &data) override;
	
	  /**
	   * @brief Returns a vector with all the variable indexes requested for post processing
	   * and an the total number of equations
	   */
	  void GetPostProcessVarIndexList(TPZVec<int> & varIndexList);
    
    /**
     * @brief Return the name of the ith postproc variable
     */
    void GetPostProcVarName(int64_t index, std::string &varname)
    {
        varname = fVars[index].fName;
        
    }

	  /**
	   * @brief Informs the vector with all the variable indexes requested for post processing
	   * and the reference to the reference postprocessed material such that the other
	   * information may be acquired.
	   */
	  void SetPostProcessVarIndexList(TPZVec<std::string> & varIndexNames, TPZMaterial * pRefMat);
	
protected:
		
		TPZManVector<TPZPostProcVar,20> fVars;
		int fDimension;
};

#endif
