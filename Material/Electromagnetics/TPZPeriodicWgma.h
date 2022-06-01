#ifndef _TPZPLANARWGMODALANALYSIS_H_
#define _TPZPLANARWGMODALANALYSIS_H_

#include "TPZScalarField.h"
#include "TPZMatGeneralisedEigenVal.h"

/**
 * @ingroup material
 * @brief This class implements the weak statement for the modal analysis of periodic planar waveguides using H1 elements.
 * It can either approximate TE or TM modes of a 2D waveguide with a 1D cross-section.
 * It also supports periodicity in the z-direction (as in photonic crystals).
 * For two-dimensional periodic waveguides, it solves the generalised
 * eigenvalue problem
   K(beta)x = beta^2 Mx
 * @note CURRENTLY ONLY IMPLEMENTED FOR DIM=2
 Formulation taken from: 
 Yasuhide Tsuji and Masanori Koshiba, "Finite Element Method Using Port Truncation by Perfectly Matched Layer Boundary Conditions for Optical Waveguide Discontinuity Problems," J. Lightwave Technol. 20, 463- (2002) 
*/
class TPZPeriodicWgma
  : public TPZScalarField,
    public TPZMatGeneralisedEigenVal {

public:

  //! Parent class's constructor
  using TPZScalarField::TPZScalarField;
  
  TPZPeriodicWgma *NewMaterial() const override;

  std::string Name() const override { return "TPZPeriodicWgma"; }

  //! Sets the propagation constant to be used in the K matrix
  void SetBeta(const CSTATE beta){fBeta = beta;}
  
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
    
  void Read(TPZStream& buf, void* context) override;

  void Write(TPZStream& buf, int withclassid) const override;


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
  void ContributeBCInternalA(const CSTATE coeffGradX,
                            const TPZMaterialDataT<CSTATE> &data, REAL weight,
                            TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef,
                            TPZBndCondT<CSTATE> &bc);
  //! Actual Contribute method (both TE/TM) for matrix B
  void ContributeInternalB(const CSTATE cGradX, const CSTATE cGradY,
                          const CSTATE cScal,
                          const TPZMaterialDataT<CSTATE> &data, REAL weight,
                          TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef);
  
  //! Propagation constant
  CSTATE fBeta;
  //! Default constructot
  TPZPeriodicWgma() = default;
};

#endif /* _TPZPLANARWGMODALANALYSIS_H_ */
