/**
 * @file TPZScattering.h
 * @brief Header file for class TPZScattering.\n
 */

#ifndef TPZSCATTERING_H
#define TPZSCATTERING_H


#include "TPZMatBase.h"
#include "TPZMatSingleSpace.h"
/**
 * @ingroup material
 * @brief This class implements the weak statement for the scattering analysis of 2D
 waveguides using HCurl elements
 * This source is implemented through the material TPZScatteringSrc.
*/
class  TPZScattering  :
  public TPZMatBase<CSTATE,TPZMatSingleSpaceT<CSTATE>>
{
  using TBase = TPZMatBase<CSTATE,TPZMatSingleSpaceT<CSTATE>>;
public:
  /**
       @brief Constructor taking a few material parameters
       @param[in] id Material identifier.
       @param[in] er Relative permittivity.
       @param[in] ur Relative permeability.
       @param[in] scale Scale for geometric domain.
       @note the `scale` param might help with floating point arithmetics on really small domains.
    */
  TPZScattering(int id, const CSTATE er, const CSTATE ur, const STATE lambda,
          const REAL scale = 1.);

  TPZScattering(int id, const TPZFMatrix<CSTATE>& er, const TPZFMatrix<CSTATE>& ur, const STATE lambda,
          const REAL scale = 1.);
  
  TPZScattering * NewMaterial() const override;
    
  std::string Name() const override { return "TPZScattering"; }

  /** @brief Returns the integrable dimension of the material */
  int Dimension() const override {return 3;}

  [[nodiscard]] int NStateVariables() const override{return 1;}


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
  
  /**
     @name ParamMethods
     @{
  */
  //! Sets the wavelength being analysed
  void SetWavelength(STATE lambda) {fLambda = lambda;}
  //! Gets the current wavelength
  [[nodiscard]] inline STATE GetWavelength() const{ return fLambda;}
  //! Sets the permeability of the material
  void SetPermeability(CSTATE ur);
  //! Sets the permeability of the material
  void SetPermeability(const TPZFMatrix<CSTATE> &ur);
  //! Gets the permeability of the material
  virtual void GetPermeability([[maybe_unused]] const TPZVec<REAL> &x,
                               TPZFMatrix<CSTATE> &ur) const;
  //! Sets the permittivity of the material
  void SetPermittivity(CSTATE er);
  //! Sets the permittivity of the material
  void SetPermittivity(const TPZFMatrix<CSTATE> &er);
  //! Gets the permittivity of the material
  virtual void GetPermittivity([[maybe_unused]] const TPZVec<REAL> &x,
                               TPZFMatrix<CSTATE> &er) const;
  /**@}*/
  
  /**
     @name ContributeMethods
     @{
  */
  void Contribute(const TPZMaterialDataT<CSTATE> &data, REAL weight,
                  TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef) override;

  void Contribute(const TPZMaterialDataT<CSTATE> &data, REAL weight,
                  TPZFMatrix<CSTATE> &ef) override {}//nothing to be done
  
  void ContributeBC(const TPZMaterialDataT<CSTATE> &data, REAL weight,
                    TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef,
                    TPZBndCondT<CSTATE> &bc) override;
  /**@}*/

  [[nodiscard]] int ClassId() const override;
protected:
  TPZScattering() = default;
  //! Relative magnetic permeability (xx, yy, zz)
  TPZFNMatrix<9,CSTATE> fUr{{1.,0,0},{0,1,0},{0,0,1}};
  //! Relative electric permittivity (xx, yy, zz)
  TPZFNMatrix<9,CSTATE> fEr{{1.,0,0},{0,1,0},{0,0,1}};
  //! Wavelength being analysed
  STATE fLambda{1.55e-9};
  //! Scale factor for the domain (helps with floating point arithmetic on small domains)
  const REAL fScaleFactor{1.};
};

#endif
