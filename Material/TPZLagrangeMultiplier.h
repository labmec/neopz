/**
 * @file TPZLagrangeMultiplier.h
 * @brief Contains the TPZLagrangeMultiplier class which implements a lagrange multiplier to be used in FEM formulations.
 */

#ifndef TPZLAGRANGEMULTIPLIER_H
#define TPZLAGRANGEMULTIPLIER_H

#include "TPZMatBase.h"
#include "TPZMatSingleSpace.h"
#include "TPZMatInterfaceSingleSpace.h"

///! Dummy class for identifying Lagrange Multiplier materials
class TPZLagrangeMultiplierBase{
    virtual bool IsLagrangeMult() = 0;
};

/// Material which implements a Lagrange Multiplier for single space materials.
template<class TVar = STATE>
class TPZLagrangeMultiplier :
    public TPZMatBase<TVar,
                      TPZMatSingleSpaceT<TVar>,
                      TPZMatInterfaceSingleSpace<TVar>>,
    public TPZLagrangeMultiplierBase
{
    using TBase = TPZMatBase<TVar,
                             TPZMatSingleSpaceT<TVar>,
                             TPZMatInterfaceSingleSpace<TVar>>;
    
    bool IsLagrangeMult() override{return true;};
    public :
	/** @brief Simple constructor */
	TPZLagrangeMultiplier() : 
        TPZRegisterClassId(&TPZLagrangeMultiplier::ClassId),
        TBase()
    {
        
    }
	/** @brief Constructor with the index of the material object within the vector */
	TPZLagrangeMultiplier(int nummat, int dimension, int nstate=1) :
        TBase(nummat),
        fNStateVariables(nstate), fDimension(dimension)
    {
        
    }
    
    TPZMaterial *NewMaterial() const override
    {
        return new TPZLagrangeMultiplier(*this);
    }
    
    /** @brief Returns the integrable dimension of the material */
    int Dimension() const override
    {return fDimension;}

    //! @name Lagrange
    /** @{*/
    
    //! Sets the multiplier.
    virtual void SetMultiplier(const TVar mult)
    {fMultiplier = mult;}
    //! Gets the multiplier.
    TVar Multiplier()
    {return fMultiplier;}
    /** @}*/
	std::string Name() const override
    {return "TPZLagrangeMultiplier";}
    
    int NStateVariables() const override
    {return fNStateVariables;}

    void FillDataRequirementsInterface(TPZMaterialData &data) const override
    { data.SetAllRequirements(false);}
    
    void Contribute(const TPZMaterialDataT<TVar> &data, REAL weight,
                    TPZFMatrix<TVar> &ek,
                    TPZFMatrix<TVar> &ef) override
    {}
    void ContributeBC(const TPZMaterialDataT<TVar> &data, REAL weight,
                      TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef,
                      TPZBndCondT<TVar> &bc) override
    {}
    void Solution(const TPZMaterialDataT<TVar> &data, int var,
                  TPZVec<TVar> &sol) override
    {}
    void ContributeInterface(const TPZMaterialDataT<TVar> &data,
                             const TPZMaterialDataT<TVar> &dataleft,
                             const TPZMaterialDataT<TVar> &dataright,
                             REAL weight, TPZFMatrix<TVar> &ek,
                             TPZFMatrix<TVar> &ef) override;
    void
    ContributeBCInterface(const TPZMaterialDataT<TVar> &data,
                          const TPZMaterialDataT<TVar> &dataleft, REAL weight,
                          TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef,
                          TPZBndCondT<TVar> &bc) override
    {}

    void SolutionInterface(const TPZMaterialDataT<TVar> &data,
                  const TPZMaterialDataT<TVar> &dataleft,
                  const TPZMaterialDataT<TVar> &dataright,
                  int var, TPZVec<TVar> &Solout) override
    {}

    void GetSolDimensions(uint64_t &u_len,
                                  uint64_t &du_row,
                                  uint64_t &du_col) const override
    {u_len=du_row=du_col=0;}
    
    void Print(std::ostream &out) const override;
    
    int ClassId() const override;
	
	void Write(TPZStream &buf, int withclassid) const override;

	void Read(TPZStream &buf, void *context) override;

protected:
    //! Number of state variables
    int fNStateVariables{0};
    
    //! Dimension associated with the material
    int fDimension;
    //! Associated multiplier
    TVar fMultiplier{1.};
};


extern template class TPZLagrangeMultiplier<STATE>;
extern template class TPZLagrangeMultiplier<CSTATE>;
#endif
