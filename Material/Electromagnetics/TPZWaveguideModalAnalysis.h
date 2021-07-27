/**

 * @file TPZWaveguideModalAnalysis.h
 * @brief Header file for class TPZWaveguideModalAnalysis.\n
 */

#ifndef TPZWAVEGUIDEMODALANALYSIS_H
#define TPZWAVEGUIDEMODALANALYSIS_H


#include "TPZMatBase.h"
#include "TPZMatCombinedSpaces.h"
#include "TPZMatGeneralisedEigenVal.h"
/**
 * @ingroup material
 * @brief This class implements the weak statement for the modal analysis of waveguides using HCurl and H1 elements.
 * It uses a 2D Hcurl space for the transversal components of the electric field and an 1D H1 space for the longitudinal component.
 * @note Formulation taken from: LEE, J.-F.; SUN, D.-K.; CENDES, Z.J. Full-wave analysis of dielectric waveguides using tangential vector finite elements.IEEE Transactions on Microwave Theory and Techniques, Institute of Electrical and Electronics Engineers (IEEE), v. 39, n. 8,p. 1262–1271, 1991
 */
class  TPZWaveguideModalAnalysis :
    public TPZMatBase<CSTATE,
                      TPZMatCombinedSpacesT<CSTATE>,
                      TPZMatGeneralisedEigenVal>
{
    using TBase = TPZMatBase<CSTATE,
                             TPZMatCombinedSpacesT<CSTATE>,
                             TPZMatGeneralisedEigenVal>;
public:

    /**
       @brief Constructor taking a few material parameters
       @param[in] id Material identifier.
       @param[in] ur Relative permeability.
       @param[in] er Relative permittivity.
       @param[in] scale Scale for geometric domain.
       @note the `scale` param might help with floating point arithmetics on really small domains.
    */
    TPZWaveguideModalAnalysis(int id, const CSTATE ur,
                              const CSTATE er, const STATE lambda,
                              const REAL &scale = 1.);
    /**
       @brief Constructor taking a few material parameters
       @param[in] id Material identifier.
       @param[in] ur Relative permeability (xx, yy, zz).
       @param[in] er Relative permittivity (xx, yy, zz).
       @param[in] scale Scale for geometric domain.
       @note the `scale` param might help with floating point arithmetics on really small domains.
    */
    TPZWaveguideModalAnalysis(int id,
                              const TPZVec<CSTATE> & ur,
                              const TPZVec<CSTATE> & er,
                              STATE lambda,
                              const REAL &scale = 1.);
    explicit TPZWaveguideModalAnalysis(int id);

    TPZWaveguideModalAnalysis * NewMaterial() const override;
    
    std::string Name() const override { return "TPZWaveguideModalAnalysis"; }
    
    /** @brief Returns the integrable dimension of the material */
    int Dimension() const override {return 2;}

    [[nodiscard]] int NStateVariables() const override{return 1;}
    
    /** @brief Index of the HCurl approximation space*/
    [[nodiscard]] inline static constexpr int HCurlIndex() { return fHCurlMeshIndex;}
    /** @brief Index of the H1 approximation space*/
    [[nodiscard]] inline static constexpr int H1Index() { return fH1MeshIndex;}

    /**
       @name ParamMethods
       @{
     */
    //! Sets the wavelength being analysed
    void SetWavelength(STATE lambda);
    //! Gets the current wavelength
    [[nodiscard]] inline STATE GetWavelength() const{ return fLambda;}
    //! Sets the permeability of the material
    void SetPermeability(CSTATE ur);
    //! Sets the permeability of the material
    void SetPermeability(const TPZVec<CSTATE> &ur);
    //! Gets the permeability of the material
    virtual void GetPermeability([[maybe_unused]] const TPZVec<REAL> &x,
                                 TPZVec<CSTATE> &ur) const;
    //! Sets the permittivity of the material
    void SetPermittivity(CSTATE ur);
    //! Sets the permittivity of the material
    void SetPermittivity(const TPZVec<CSTATE> &er);
    //! Gets the permittivity of the material
    virtual void GetPermittivity([[maybe_unused]] const TPZVec<REAL> &x,
                                 TPZVec<CSTATE> &er) const;
    /**@}*/
    /**
       @name ContributeMethods
       @{
    */
    void Contribute(const TPZVec<TPZMaterialDataT<CSTATE>> &datavec, REAL weight,
                    TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef) override;
    
    void ContributeBC(const TPZVec<TPZMaterialDataT<CSTATE>> &datavec, REAL weight,
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
    void Solution(const TPZVec<TPZMaterialDataT<CSTATE>> &datavec,
                  int var, TPZVec<CSTATE> &solout) override;
    //! Get the propagation constant used for post processing the solution
    inline const CSTATE &GetKz() const
    {return fKz;}
    //! Set the propagation constant used for post processing the solution
    inline void SetKz(const CSTATE &kz)
    { fKz = kz;}
    /** @brief Whether to print only the real part of the electromagnetic field.
     If false, it prints the magnitude of the field.*/
    [[nodiscard]] inline bool ShouldPrintFieldRealPart() const
    { return fPrintFieldRealPart;}
    /** @brief Set to print only the real part of the electromagnetic field.
     If set to false, it prints the magnitude of the field*/
    void SetPrintFieldRealPart(bool printFieldRealPart)
    { fPrintFieldRealPart = printFieldRealPart;}
    /**@}*/

    /** @name GeneralisedMethods */
    /** @{*/
    //! Set the material to compute matrix A
    void SetMatrixA() override;
    //! Set the material to compute matrix B
    void SetMatrixB() override;
    /**@{*/
protected:
    /** @name InternalGeneralisedTypes */
    /** @{ */
    using TContributeType =
        std::function<void (const TPZVec<TPZMaterialDataT<CSTATE>> &datavec,
                            REAL weight,
                            TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef,
                            const TPZVec<CSTATE> & er, const TPZVec<CSTATE> & ur)>;
    using TContributeBCType =
        std::function<void (const TPZVec<TPZMaterialDataT<CSTATE>> &datavec,
                            REAL weight,
                            TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef,
                            TPZBndCondT<CSTATE> &bc)>;
    /** @} */
    //! Relative magnetic permeability (xx, yy, zz)
    TPZManVector<CSTATE,3> fUr{{1.,1.,1.}};
    //! Relative electric permittivity in the x direction
    TPZManVector<CSTATE,3> fEr{{1.,1.,1.}};
    //! Wavelength being analysed
    STATE fLambda{1.55e-9};
    //! Scale factor for the domain (helps with floating point arithmetic on small domains)
    const REAL fScaleFactor{1.};
    //! Alternates between printing the real part of magnitude of the field
    bool fPrintFieldRealPart{true};
    //!Fixes a propagation constant for printing the solution
    CSTATE fKz{1};
    static constexpr int fH1MeshIndex{1};
    static constexpr int fHCurlMeshIndex{0};

    /** @name InternalGeneralisedAttributes */
    /** @{ */
    //! Pointer to the current Contribute function
    TContributeType fCurrentContribute{nullptr};
    //! Pointer to the current ContributeBC function
    TContributeBCType fCurrentContributeBC{nullptr};
    /** @} */
    TPZWaveguideModalAnalysis();

    /** @name InternalContributeMethods */
    /** @{*/
    //! Contribution of the A Matrix
    void ContributeA(const TPZVec<TPZMaterialDataT<CSTATE>> &datavec, REAL weight,
                     TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef,
                     const TPZVec<CSTATE> & er, const TPZVec<CSTATE> & ur);
    //! Boundary contribution of the A Matrix
    void ContributeBCA(const TPZVec<TPZMaterialDataT<CSTATE>> &datavec, REAL weight,
                      TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef,
                      TPZBndCondT<CSTATE> &bc);
    //! Contribution of the B Matrix
    void ContributeB(const TPZVec<TPZMaterialDataT<CSTATE>> &datavec, REAL weight,
                     TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef,
                     const TPZVec<CSTATE> & er, const TPZVec<CSTATE> & ur);
    //! Boundary contribution of the B Matrix
    void ContributeBCB(const TPZVec<TPZMaterialDataT<CSTATE>> &datavec, REAL weight,
                      TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef,
                      TPZBndCondT<CSTATE> &bc);
    /** @{ */
};

#endif
