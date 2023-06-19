/**

 * @file TPZWgma.h
 * @brief Header file for class TPZWgma.\n
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
class  TPZWgma :
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
       @param[in] er Relative permittivity.
       @param[in] ur Relative permeability.
       @param[in] scale Scale for geometric domain.
       @note the `scale` param might help with floating point arithmetics on really small domains.
    */
    TPZWgma(int id, const CSTATE er,
            const CSTATE ur, const STATE lambda,
            const REAL &scale = 1.);
    /**
       @brief Constructor taking a few material parameters
       @param[in] id Material identifier.
       @param[in] er Relative permittivity (xx, yy, zz).
       @param[in] ur Relative permeability (xx, yy, zz).
       @param[in] scale Scale for geometric domain.
       @note the `scale` param might help with floating point arithmetics on really small domains.
    */
    TPZWgma(int id,
            const TPZFMatrix<CSTATE> & er,
            const TPZFMatrix<CSTATE> & ur,
            STATE lambda,
            const REAL &scale = 1.);
    explicit TPZWgma(int id);

    TPZWgma * NewMaterial() const override;
    
    std::string Name() const override { return "TPZWgma"; }
    
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
    void SetPermeability(const TPZFMatrix<CSTATE> &ur);
    //! Gets the permeability of the material
    virtual void GetPermeability([[maybe_unused]] const TPZVec<REAL> &x,
                                 TPZFMatrix<CSTATE> &ur) const;
    //! Sets the permittivity of the material
    void SetPermittivity(CSTATE ur);
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
    void Contribute(const TPZVec<TPZMaterialDataT<CSTATE>> &datavec, REAL weight,
                    TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef) override;
    
    void ContributeBC(const TPZVec<TPZMaterialDataT<CSTATE>> &datavec, REAL weight,
                      TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef,
                      TPZBndCondT<CSTATE> &bc) override;
    /**@}*/
    /**
       @name SolutionMethods
       @{*/
    /** @brief Variable index of a given solution.
        Possibilities are:
        -Et_real
        -Ez_real
        -Et_abs
        -Ez_abs
        -Material
        -P_H1
        -P_HCurl
    */
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
                            const TPZFMatrix<CSTATE> & er,
                            const TPZFMatrix<CSTATE> & ur)>;
    using TContributeBCType =
        std::function<void (const TPZVec<TPZMaterialDataT<CSTATE>> &datavec,
                            REAL weight,
                            TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef,
                            TPZBndCondT<CSTATE> &bc)>;
    /** @} */
    //! Relative magnetic permeability
    TPZFNMatrix<9,CSTATE> fUr{{1.,0,0},{0,1,0},{0,0,1}};
    //! Relative electric permittivity
    TPZFNMatrix<9,CSTATE> fEr{{1.,0,0},{0,1,0},{0,0,1}};
    //! Wavelength being analysed
    STATE fLambda{1.55e-9};
    //! Scale factor for the domain (helps with floating point arithmetic on small domains)
    const REAL fScaleFactor{1.};
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
    TPZWgma();//< Default constructor

    /** @name InternalContributeMethods */
    /** @{*/
    //! Contribution of the A Matrix
    void ContributeA(const TPZVec<TPZMaterialDataT<CSTATE>> &datavec, REAL weight,
                     TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef,
                     const TPZFMatrix<CSTATE> & er, const TPZFMatrix<CSTATE> & ur);
    //! Boundary contribution of the A Matrix
    void ContributeBCA(const TPZVec<TPZMaterialDataT<CSTATE>> &datavec, REAL weight,
                      TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef,
                      TPZBndCondT<CSTATE> &bc);
    //! Contribution of the B Matrix
    void ContributeB(const TPZVec<TPZMaterialDataT<CSTATE>> &datavec, REAL weight,
                     TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef,
                     const TPZFMatrix<CSTATE> & er, const TPZFMatrix<CSTATE> & ur);
    //! Boundary contribution of the B Matrix
    void ContributeBCB(const TPZVec<TPZMaterialDataT<CSTATE>> &datavec, REAL weight,
                      TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef,
                      TPZBndCondT<CSTATE> &bc);
    /** @{ */
};

#endif
