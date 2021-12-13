/**
 * @file TPZMatPoisson.h
 * @brief Contains the TPZMatPoisson class which a H1 formulation of the Poisson equation
 */

#ifndef TPZMATPOISSON_H
#define TPZMATPOISSON_H

#include "TPZMatBase.h"
#include "TPZMatSingleSpace.h"
#include "TPZMatLoadCases.h"
#include "TPZMatErrorSingleSpace.h"


/**
 * @brief Implements a H1 formulation of the poisson equation.
 * This material uses a scalar approximation space in a geometric domain of
 * dimension defined by the user (1, 2 or 3D).
 * It can solve for multiple values of rhs at once, but for calculating the error,
 * it will only consider one solution at a time.
 */
template<class TVar=STATE>
class TPZMatPoisson :
    public TPZMatBase<TVar,
                      TPZMatSingleSpaceT<TVar>,
                      TPZMatErrorSingleSpace<TVar>,
                      TPZMatLoadCases<TVar>>{
	using TBase = TPZMatBase<TVar,
                             TPZMatSingleSpaceT<TVar>,
                             TPZMatErrorSingleSpace<TVar>,
                             TPZMatLoadCases<TVar>>;
public:
    //! Default constructor
    TPZMatPoisson() = default;
    /**
	 * @brief Class constructor 
	 * @param id material id
	 * @param dim problem dimension
	 * @param nstate number of state variables
	 */
	TPZMatPoisson(int id, int dim);

    std::string Name() const override { return "TPZMatPoisson"; }
	
	/** @brief Solution indices of post-processing */
	enum ESolutionVars { ENone = 0, ESolution = 1 , EDerivative = 2};
	
    /**
     * @brief Set a scale factor for the stiffness matrix and right hand side
     * the default value of the scale factor is 1
     */
    void SetScaleFactor(REAL scale)
    {fScale = scale;}
    
    [[nodiscard]] REAL ScaleFactor() const
    {return fScale;}

	int Dimension() const  override { return this->fDim; }

    int NStateVariables() const override{ return 1;}
    //H1 norm L2 norm H1 seminorm
    int NEvalErrors() const override { return 3;}
    /** @brief Sets problem dimension */
    virtual void SetDimension(int dim) { this->fDim = dim; }

    void Contribute(const TPZMaterialDataT<TVar> &data,
                    REAL weight,
                    TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef) override;

    void ContributeBC(const TPZMaterialDataT<TVar> &data, REAL weight,
                      TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef,
                      TPZBndCondT<TVar> &bc) override;
    /** @brief To create another material of the same type */
	TPZMaterial * NewMaterial() const override;
	
	/** @brief It returns the variable index associated with the name */
	int VariableIndex(const std::string &name)const override;
	
	int NSolutionVariables(int var) const override;

    void Solution(const TPZMaterialDataT<TVar> &data,
                  int var, TPZVec<TVar> &solOut) override;
    
    void GetSolDimensions(uint64_t &u_len,
                          uint64_t &du_row,
                          uint64_t &du_col) const override;

    void Errors(const TPZMaterialDataT<TVar>&data,
                TPZVec<REAL> &errors) override;
    
    virtual int ClassId() const override;
protected:
	
	/** @brief Problem dimension */
	int fDim;
	
	/** @brief Constant solution vector. Ignored if forcing function is set. */
	TPZVec<TVar> fSol;
    
    /** @brief Scale factor applied to the stiffness matrix and right hand side */
    REAL fScale{1};
};


extern template class TPZMatPoisson<STATE>;
extern template class TPZMatPoisson<CSTATE>;
#endif
