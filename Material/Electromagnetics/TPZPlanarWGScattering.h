/**
 * @file TPZPlanarWGScattering.h
 * @brief Header file for class TPZPlanarWGScattering.\n
 */

#ifndef TPZPLANARSCATTERING_H
#define TPZPLANARSCATTERING_H


#include "TPZMatBase.h"
#include "TPZMatSingleSpace.h"
/**
 * @ingroup material
 * @brief This class implements the weak statement for the scattering analysis of planar waveguides using H1 elements.
 * It can either approximate TE or TM modes of a 2D waveguide with a 1D cross-section.
 * Currently, only sources aligned with the y direction of a 2D domain can be used.
 * This source is implemented through a ForcingFunctionBCType. The matrix-type
 * parameter is not used, and the vector-type is expected to contain source-value
 * in the first position and beta-value in the second position.
 * @note Formulation taken from: 
 Yasuhide Tsuji and Masanori Koshiba, "Finite Element Method Using Port Truncation by Perfectly Matched Layer Boundary Conditions for Optical Waveguide Discontinuity Problems," J. Lightwave Technol. 20, 463- (2002) 
*/
class  TPZPlanarWGScattering :
  public TPZMatBase<CSTATE,TPZMatSingleSpaceT<CSTATE>>
{
  using TBase = TPZMatBase<CSTATE,TPZMatSingleSpaceT<CSTATE>>;
public:
  //! Sets for which mode type the problem will be solved.
  enum class ModeType{TE, TM};
  /**
     @brief Constructor taking a few material parameters
     @param[in] id Material identifier.
     @param[in] ur Relative permeability.
     @param[in] er Relative permittivity.
     @param[in] scale Scale for geometric domain.
     @note the `scale` param might help with floating point arithmetics on really small domains.
  */
  TPZPlanarWGScattering(int id, const CSTATE ur,const CSTATE er,
                        const STATE lambda,
                        const ModeType mode,
                        const REAL &scale = 1.);

  TPZPlanarWGScattering * NewMaterial() const override;
    
  std::string Name() const override { return "TPZPlanarWGScattering"; }
  
  /** @brief Returns the integrable dimension of the material */
  int Dimension() const override {return 2;}

  [[nodiscard]] int NStateVariables() const override{return 1;}  

  /**
     @name ParamMethods
     @{
  */
  //! Sets the wavelength being analysed
  void SetWavelength(STATE lambda);
  //! Gets the current wavelength
  [[nodiscard]] inline STATE GetWavelength() const{ return fLambda;}
  //! Sets the permeability of the material
  void SetPermeability(const CSTATE ur);
  //! Gets the permeability of the material
  virtual TPZVec<CSTATE> GetPermeability([[maybe_unused]] const TPZVec<REAL> &x) const;
  //! Sets the permittivity of the material
  void SetPermittivity(const CSTATE er);
  //! Gets the permittivity of the material
  virtual TPZVec<CSTATE> GetPermittivity([[maybe_unused]] const TPZVec<REAL> &x) const;
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
  /**
     @name SolutionMethods
     @{*/
  //! Variable index of a given solution
  int VariableIndex(const std::string &name) const override;
  //! Number of variables associated with a given solution
  int NSolutionVariables(int var) const override;
  //! Computes the solution at an integration point
  void Solution(const TPZMaterialDataT<CSTATE> &data,
                int var, TPZVec<CSTATE> &solout) override;
  /**@}*/
protected:
  //! Actual Contribute method (both TE/TM)
  void ContributeInternal(const CSTATE cGradX, const CSTATE cGradY, const CSTATE cScal,
                          const TPZMaterialDataT<CSTATE> &data, REAL weight,
                          TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef);
  //! Actual ContributeBC method (both TE/TM)
  void ContributeBCInternal(const CSTATE coeffGradX,
                    const TPZMaterialDataT<CSTATE> &data, REAL weight,
                    TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef,
                    TPZBndCondT<CSTATE> &bc);
  //! Mode type being solved
  ModeType fMode{ModeType::TE};
  //! Relative magnetic permeability
  CSTATE fUr{1};
  //! Relative electric permittivity
  CSTATE fEr{1};
  //! Wavelength being analysed
  STATE fLambda{1.55e-9};
  //! Scale factor for the domain (helps with floating point arithmetic on small domains)
  const REAL fScaleFactor{1.};

  TPZPlanarWGScattering();
};

#endif
