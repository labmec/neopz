/**
 * @file
 * @brief Contains the TPZMatLaplacian class.
 */

#ifndef MATLAPLACLAGRANGEDH
#define MATLAPLACLAGRANGEDH

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
class TPZMatLaplacianLagrange : public TPZMatLaplacian {
	
	protected :
	
public:
	

	TPZMatLaplacianLagrange(int nummat, int dim);

  TPZMatLaplacianLagrange(int matid) : TPZMatLaplacian(matid)
  {

  }

	TPZMatLaplacianLagrange();

	TPZMatLaplacianLagrange(const TPZMatLaplacianLagrange &copy) : TPZMatLaplacian(copy)
    {
        
    }

	virtual ~TPZMatLaplacianLagrange();

	TPZMatLaplacianLagrange &operator=(const TPZMatLaplacianLagrange &copy);


	virtual TPZMaterial * NewMaterial(){
		return new TPZMatLaplacianLagrange(*this);
	}


	virtual void Print(std::ostream & out);

	virtual std::string Name() { return "TPZMatLaplacianLagrange"; }

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

    void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout);
    
    virtual int VariableIndex(const std::string &name);
    virtual int NSolutionVariables(int var);


  public:



    virtual int ClassId() const {
        return TPZMatLaplacianLagrangeID;
    }

	virtual void Write(TPZStream &buf, int withclassid);

	virtual void Read(TPZStream &buf, void *context);

};

#endif

