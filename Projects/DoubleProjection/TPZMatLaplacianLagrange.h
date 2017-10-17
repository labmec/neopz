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

/**
 * @ingroup  Material class of the Discontinuous Petrov-Galerkin (DPG) method
 * @ingroup  Find: (u,q) in Xh and e in Vh, so that
 * @brief \f$ (e,v) + b((u,q),v) = f(v), for all v in Vh \f$
 * @brief  \f$        b((z,r),e) = 0, for all (z,r) in Xh \f$
 */
class TPZMatLaplacianLagrange : public TPZMatLaplacian {
	
	protected :
    
    //flag used to change operator
    /*
     * Find: (u,q) in Xh and e in Vh, so that
     * b(e,v) + ((u,q),v) = f(v), for all v in Vh
     *          ((z,r),e) = 0,    for all (z,r) in Xh
     */
    bool fIsDPGPhil;
	
    //use the douple projection method (MDP)
    /* Vr: rich space
     * Find: u in Xh and ur in Vr, so that
     * b(ur,vr) = f(vr), for all vr in Vr
     * (ur,v) - (u,v) = 0,    for all v in Xh
     */
    bool fUseMDP;
    
public:
	

  TPZMatLaplacianLagrange(int nummat, int dim);

  TPZMatLaplacianLagrange(int matid) : TPZMatLaplacian(matid)
  {

  }

	TPZMatLaplacianLagrange();

	TPZMatLaplacianLagrange(const TPZMatLaplacianLagrange &copy);
    
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
     * @brief Contribute to DPG method
     * @param datavec [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     */
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec,REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
    //Contribute to DPG method by Philippe
    void ContributeDPGPhil(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    //contibute to Douple Projection Method
    void ContributeMDP(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    /** @} */

    virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout);

    virtual int VariableIndex(const std::string &name);
    
	virtual int NSolutionVariables(int var);
    
    void SetDPGPhil(){
        fIsDPGPhil = true;
    }
    
    void SetMDP(){
        fUseMDP = true;
    }
    
  public:



    virtual int ClassId() const{
        return TPZMatLaplacianLagrangeID;
    }

	virtual void Write(TPZStream &buf, int withclassid) const;

	virtual void Read(TPZStream &buf, void *context);

};

#endif

