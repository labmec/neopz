#ifndef _TPZPLANARWGMA_H_
#define _TPZPLANARWGMA_H_

#include "TPZScalarField.h"
#include "TPZMatGeneralisedEigenVal.h"

/**
 * @ingroup material
 * @brief This class implements the weak statement for the modal analysis of planar waveguides using H1 elements.
 * It can either approximate TE or TM modes of a 2D waveguide with a 1D cross-section.
 * It also supports periodicity in the z-direction (as in photonic crystals).
 * For one-dimensional cross sections, it solves the generalised
 * eigenvalue problem
   Kx = beta^2 Mx
 * @note CURRENTLY ONLY IMPLEMENTED FOR DIM=1
*/
class TPZPlanarWgma
  : public TPZScalarField,
    public TPZMatGeneralisedEigenVal {

public:

  //! Parent class's constructor
  using TPZScalarField::TPZScalarField;
  
  TPZPlanarWgma *NewMaterial() const override;

  std::string Name() const override { return "TPZPlanarWgma"; }


  /** @brief Returns the integrable dimension of the material */
  int Dimension() const override {return 1;}
  
  /**@}*/
  /**
     @name ContributeMethods
     @{
  */
  void Contribute(const TPZMaterialDataT<CSTATE> &data, REAL weight,
                  TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef) override;

  void ContributeBC(const TPZMaterialDataT<CSTATE> &data, REAL weight,
                    TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef,
                    TPZBndCondT<CSTATE> &bc) override;
  /**@}*/

  [[nodiscard]] int ClassId() const override;
  
  void Read(TPZStream& buf, void* context) override{
    TPZScalarField::Read(buf,context);
  }

  void Write(TPZStream& buf, int withclassid) const override{
    TPZScalarField::Write(buf,withclassid);
  }
  
  TPZBndCondT<CSTATE>* CreateBC(TPZMaterial *reference,
                                int id, int type,
                                const TPZFMatrix<CSTATE> &val1,
                                const TPZVec<CSTATE> &val2) override
  {
    return new  TPZBndCondBase<CSTATE,
                               TPZMatSingleSpaceBC<CSTATE>,
                               TPZMatGeneralisedEigenValBC>
      (reference,id, type,val1,val2);
  }
  
protected:
  //! Actual Contribute method (both TE/TM) for matrix A
  void ContributeInternalA(const CSTATE cGradX, const CSTATE cGradY,
                           const CSTATE cScal,
                           const TPZMaterialDataT<CSTATE> &data, REAL weight,
                           TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef);
  //! Actual ContributeBC method (both TE/TM) for matrix A
  void ContributeBCInternalA(const TPZMaterialDataT<CSTATE> &data, REAL weight,
                             TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef,
                             TPZBndCondT<CSTATE> &bc);
  //! Actual Contribute method (both TE/TM) for matrix B
  void ContributeInternalB(const CSTATE cGradX, const CSTATE cGradY,
                           const CSTATE cScal,
                           const TPZMaterialDataT<CSTATE> &data, REAL weight,
                           TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef);
  //! Default constructot
  TPZPlanarWgma() = default;
};

#endif /* _TPZPLANARWGMA_H_ */
