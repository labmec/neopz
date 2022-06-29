#ifndef _TPZSCALARFIELD_H_
#define _TPZSCALARFIELD_H_
#include "TPZMatBase.h"
#include "TPZMatSingleSpace.h"

class  TPZScalarField :
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
  TPZScalarField(int id, const CSTATE er,const CSTATE ur,
                 const STATE lambda,
                 const ModeType mode,
                 const REAL scale = 1.);

  std::string Name() const override { return "TPZScalarField"; }
  
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
  virtual void GetPermeability([[maybe_unused]] const TPZVec<REAL> &x,
                               TPZVec<CSTATE> &ur) const;
  //! Sets the permittivity of the material
  void SetPermittivity(const CSTATE er);
  //! Gets the permittivity of the material
  virtual void GetPermittivity([[maybe_unused]] const TPZVec<REAL> &x,
                               TPZVec<CSTATE> &er) const;
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


  [[nodiscard]] int ClassId() const override;
    
  void Read(TPZStream& buf, void* context) override;

  void Write(TPZStream& buf, int withclassid) const override;
  
protected:
  //! Mode type being solved
  ModeType fMode{ModeType::TE};
  //! Relative magnetic permeability
  CSTATE fUr{1};
  //! Relative electric permittivity
  CSTATE fEr{1};
  //! Wavelength being analysed
  STATE fLambda{1.55e-9};
  //! Scale factor for the domain (helps with floating point arithmetic on small domains)
  REAL fScaleFactor{1.};

  TPZScalarField() = default;
};

#endif /* _TPZSCALARFIELD_H_ */
