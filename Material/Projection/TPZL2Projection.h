/**
 * @file TPZL2Projection.h
 * @brief Contains the TPZL2Projection class which implements an L2 projection of a given scalar solution
 */

#ifndef TPZL2PROJECTION_H
#define TPZL2PROJECTION_H

#include "TPZMatBase.h"
#include "TPZMatSingleSpace.h"
#include "TPZMatErrorSingleSpace.h"


/**
 * @brief Implements an L2 projection of a given solution on a scalar approximation space.
 * This material uses a scalar approximation space in a geometric domain of
 * dimension defined by the user (1, 2 or 3D).
 * It can project several solutions at once, but for calculating the error,
 * it will only consider one solution at a time.
 * The solutions to be projected are given by the Forcing Function of the 
 * TPZMaterialT.  It is expected to be a TPZVec<TVar> of size `nsol`.
 * For the error analysis, the solution is set through the ExactSol of
 * TPZMatErrorSingleSpace.
 */
template<class TVar=STATE>
class TPZL2Projection :
    public TPZMatBase<TVar,TPZMatSingleSpaceT<TVar>,TPZMatErrorSingleSpace<TVar>>{
	using TBase = TPZMatBase<TVar,TPZMatSingleSpaceT<TVar>,TPZMatErrorSingleSpace<TVar>>;	
public:
    //! Default constructor
    TPZL2Projection() = default;
    /**
	 * @brief Class constructor 
	 * @param id material id
	 * @param dim problem dimension
	 * @param nstate number of state variables
	 */
	TPZL2Projection(int id, int dim, int nstate=1);
	/**
	 * @brief Class constructor setting a default constant solution
	 * @param id material id
	 * @param dim problem dimension
	 * @param nstate number of state variables
	 * @param sol constant solution. Ignored if a forcing function is set.
	 */
	TPZL2Projection(int id, int dim, int nstate, const TPZVec<TVar> &sol);

    std::string Name() const override { return "TPZL2Projection"; }
	
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
    
    /** @brief Sets problem dimension */
    virtual void SetDimension(int dim) { this->fDim = dim; }
	
    int NStateVariables() const override { return this->fNStateVars; }
    void Contribute(const TPZMaterialDataT<TVar> &data,
                    REAL weight,
                    TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef) override;

    void ContributeBC(const TPZMaterialDataT<TVar> &data, REAL weight,
                      TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef,
                      TPZBndCondT<TVar> &bc) override;
    /** @brief To create another material of the same type */
	TPZMaterial * NewMaterial() const override;
	
	/** @brief It returns the variable index associated with the name */
	int VariableIndex(const std::string &name) const override;
	
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
	
	/** @brief Number of state variables */
	int fNStateVars{1};
	
	/** @brief Constant solution vector. Ignored if forcing function is set. */
	TPZVec<TVar> fSol;
    
    /** @brief Scale factor applied to the stiffness matrix and right hand side */
    REAL fScale{1};
};


extern template class TPZL2Projection<STATE>;
extern template class TPZL2Projection<CSTATE>;
#endif
