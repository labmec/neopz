/**
 * @file
 * @brief Contains the TPZMatLaplacian class.
 */

#ifndef MATLAPLACHYBRIDH
#define MATLAPLACHYBRIDH

#include <iostream>
#include "pzdiscgal.h"
#include "pzfmatrix.h"
#include "TPZMatLaplacian.h"

/**
 * @ingroup material
 * @brief \f$ -fK Laplac(u) = fXf  \f$
 */
/**
 * \f$ -fK Laplac(u) = fXf  \f$
 */
class TPZMatLaplacianHybrid : public TPZMatLaplacian {
	
	protected :
	
public:
	

	TPZMatLaplacianHybrid(int nummat, int dim);

  TPZMatLaplacianHybrid(int matid) : TPZMatLaplacian(matid)
  {

  }

	TPZMatLaplacianHybrid();

	TPZMatLaplacianHybrid(const TPZMatLaplacianHybrid &copy) : TPZMatLaplacian(copy)
    {
        
    }

	virtual ~TPZMatLaplacianHybrid();

	TPZMatLaplacianHybrid &operator=(const TPZMatLaplacianHybrid &copy);


	virtual TPZMaterial * NewMaterial(){
		return new TPZMatLaplacianHybrid(*this);
	}


	virtual void Print(std::ostream & out);

	virtual std::string Name() { return "TPZMatLaplacianHybrid"; }

	/**
	 * @name Contribute methods (weak formulation)
	 * @{
	 */


    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation.
     * @param datavec [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     */
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);

    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
    {
        ContributeBC(datavec[0],weight,ek,ef,bc);
    }
    
    virtual void ContributeBC(TPZMaterialData &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    

    void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout);
    
    virtual int VariableIndex(const std::string &name);
    virtual int NSolutionVariables(int var);

    

  public:

    virtual int NEvalErrors() {return 4;}


    public:
virtual int ClassId() const;


	virtual void Write(TPZStream &buf, int withclassid) const;

	virtual void Read(TPZStream &buf, void *context);

};

#endif

