/**
 * Voigt notation ordering used by TPZLinearElasticityConstitutive
 *
 * This class builds the constitutive matrix C and the strain–displacement matrix B
 * assuming the following component orderings in Voigt notation.
 *
 * 3D elasticity (dimension == 3):
 *  - Strain vector ε (6x1) and stress vector σ (6x1):
 *      [ εxx, εyy, εzz, γyz, γxz, γxy ]^T
 *      [ σxx, σyy, σzz, τyz, τxz, τxy ]^T
 *
 * 2D plane stress (dimension == 2 and fPlaneStress == true):
 *  - Strain vector ε (3x1) and stress vector σ (3x1):
 *      [ εxx, εyy, γxy ]^T
 *      [ σxx, σyy, τxy ]^T
 *
 * 2D plane strain (dimension == 2 and fPlaneStress == false):
 *  - In-plane Voigt form for assembled C and B (3x1):
 *      [ εxx, εyy, γxy ]^T
 *      [ σxx, σyy, τxy ]^T
 *    Note: εzz is constrained by plane strain (εzz = 0), which is accounted for
 *    internally in the constitutive matrix construction; out-of-plane shear
 *    components γxz and γyz are zero in pure 2D formulations.
 *
 * Shear convention:
 *  - Engineering shear strains γij = 2*εij (i ≠ j) are used in Voigt form.
 *  - Shear stresses τij are the corresponding stress components.
 *
 * Mapping to matrices:
 *  - C maps Voigt strains to Voigt stresses: σ = C * ε
 *  - B maps nodal displacements to Voigt strains: ε = B * u
 *
 * Orthotropic materials:
 *  - Material axes (principal directions) are defined by fPrincipalAxes and are
 *    rotated to the global frame via the internal rotation matrix. The Voigt
 *    ordering above applies in the global frame after rotation.
 */
#ifndef TPZLINEARELASTICITYCONSTITUTIVE_H
#define TPZLINEARELASTICITYCONSTITUTIVE_H

#include "pzreal.h"
#include "pzfmatrix.h"
#include "pzmanvector.h"

/**
 * @brief Interface class to define the constitutive model for linear elasticity
 * Provides C (constitutive) and B (strain-displacement) matrices
 * Supports orthotropic properties with principal axis directions
 */

class TPZLinearElasticityConstitutive {
    
public:
    
    /// Default constructor
    TPZLinearElasticityConstitutive();
    
    /// Destructor
    virtual ~TPZLinearElasticityConstitutive();
    
    /**
     * @brief Constructor for isotropic material
     * @param E Young's modulus
     * @param nu Poisson's ratio
     */
    TPZLinearElasticityConstitutive(REAL E, REAL nu);
    
    /**
     * @brief Constructor for orthotropic material
     * @param Ex Young's modulus in x direction
     * @param Ey Young's modulus in y direction
     * @param Ez Young's modulus in z direction
     * @param nuxy Poisson's ratio xy
     * @param nuyz Poisson's ratio yz
     * @param nuxz Poisson's ratio xz
     * @param Gxy Shear modulus xy
     * @param Gyz Shear modulus yz
     * @param Gxz Shear modulus xz
     */
    TPZLinearElasticityConstitutive(REAL Ex, REAL Ey, REAL Ez,
                                   REAL nuxy, REAL nuyz, REAL nuxz,
                                   REAL Gxy, REAL Gyz, REAL Gxz);
    
    /// @brief Constructor for orthotropic material in 2D
    /// @param Ex Young's modulus in x direction
    /// @param Ey Young's modulus in y direction
    /// @param nuxy Poisson's ratio xy
    /// @param Gxy Shear modulus xy
    /// @param planeStress If true, plane stress; if false, plane strain
    TPZLinearElasticityConstitutive(REAL Ex, REAL Ey,
                                   REAL nuxy,
                                   REAL Gxy, bool planeStress);
    /// Copy constructor
    TPZLinearElasticityConstitutive(const TPZLinearElasticityConstitutive &copy);
    
    /// Assignment operator
    TPZLinearElasticityConstitutive &operator=(const TPZLinearElasticityConstitutive &copy);
    
    /**
     * @brief Set isotropic material properties
     * @param E Young's modulus
     * @param nu Poisson's ratio
     */
    void SetIsotropicProperties(REAL E, REAL nu);
    
    /**
     * @brief Set orthotropic material properties
     */
    void SetOrthotropicProperties(REAL Ex, REAL Ey, REAL Ez,
                                  REAL nuxy, REAL nuyz, REAL nuxz,
                                  REAL Gxy, REAL Gyz, REAL Gxz);
    
    /**
     * @brief Set orthotropic material properties
     */
    void SetOrthotropicProperties(REAL Ex, REAL Ey,
                                  REAL nuxy,
                                  REAL Gxy);
    



    /**
     * @brief Convert a symmetric tensor (matrix form) to Voigt vector representation.
     * @param T Input tensor:
     *        - If fDimension==3: 3x3 symmetric matrix
     *        - If fDimension==2: 2x2 symmetric matrix (in-plane components)
     * @param v Output Voigt vector:
     *        - If fDimension==3: 6x1 [ xx, yy, zz, yz, xz, xy ]^T
     *        - If fDimension==2: 3x1 [ xx, yy, xy ]^T
     * @param engineeringShear If true, uses engineering shear convention (γij = 2*εij).
     *        Set true for strains, false for stresses.
     *
     * Notes:
     *  - For strains, set engineeringShear=true to apply factor 2 on off-diagonal terms.
     *  - For stresses, set engineeringShear=false (no factor 2).
     */
    void TensorToVoigt(const TPZFMatrix<REAL> &T, TPZFMatrix<REAL> &v, bool engineeringShear) const;

    /**
     * @brief Convert a Voigt vector to symmetric tensor (matrix form).
     * @param v Input Voigt vector:
     *        - If fDimension==3: 6x1 [ xx, yy, zz, yz, xz, xy ]^T
     *        - If fDimension==2: 3x1 [ xx, yy, xy ]^T
     * @param T Output tensor:
     *        - If fDimension==3: 3x3 symmetric matrix
     *        - If fDimension==2: 2x2 symmetric matrix (in-plane components)
     * @param engineeringShear If true, expects engineering shear convention in v (γij).
     *        Set true for strains, false for stresses.
     *
     * Notes:
     *  - For strains (engineeringShear=true), off-diagonal tensor entries are v_shear/2.
     *  - For stresses (engineeringShear=false), off-diagonal tensor entries are v_shear.
     *  - In 2D, out-of-plane components are not filled.
     */
    void VoigtToTensor(const TPZFMatrix<REAL> &v, TPZFMatrix<REAL> &T, bool engineeringShear) const;



    /**
     * @brief Compute stiffness matrix D
     * @param C Output constitutive matrix
     *
     * If dimension==2, behavior depends on fPlaneStress flag:
     *  - fPlaneStress == true  -> plane stress
     *  - fPlaneStress == false -> plane strain (plane elasticity)
     */
    virtual void ComputeStiffnessMatrix(TPZFMatrix<REAL> &D) const;
    
    /// @brief Compute the constitutive matrix for three dimensional elasticity
    /// @param D constitutive matrix
    virtual void ComputeStiffnessMatrix3D(TPZFMatrix<REAL> &D) const;

    /**
     * @brief Compute constitutive matrix for plane stress (2D)
     * @param D Output constitutive matrix
     */
    virtual void ComputeStiffnessMatrixPlaneStress(TPZFMatrix<REAL> &D) const;

    /**
     * @brief Compute constitutive matrix for plane strain (plane elasticity, 2D)
     * @param D Output constitutive matrix
     */
    virtual void ComputeStiffnessMatrixPlaneStrain(TPZFMatrix<REAL> &D) const;

    /**
     * @brief Compute constitutive matrix D for 2D and 3D in principal material directions
     * @param D Output constitutive matrix
     *
     * If dimension==2, behavior depends on fPlaneStress flag:
     *  - fPlaneStress == true  -> plane stress
     *  - fPlaneStress == false -> plane strain
     */
    virtual void ComputeStiffnessMatrix3DPrincipal(TPZFMatrix<REAL> &D) const;

    /**
     * @brief Build the 3D compliance matrix in the principal material coordinate system.
     *
     * This method assembles the 6x6 compliance matrix S (also called the flexibility
     * matrix) in Voigt notation, aligned with the principal axes defined by the material.
     * The resulting matrix relates stress to strain via:
     *    {epsilon_principal} = S_principal * {sigma_principal}
     *
     * Characteristics:
     *  - Voigt ordering is [ xx, yy, zz, yz, xz, xy ].
     *  - Matrix size is 6x6 for 3D problems.
     *  - Uses the material's orthotropic properties (Ex, Ey, Ez, nuxy, nuyz, nuxz, Gxy, Gyz, Gxz).
     *    If the material is isotropic, the orthotropic formulation degenerates to the isotropic case.
     *  - Principal directions are defined by fPrincipalAxes; this routine builds S in that local frame.
     *    No rotation to global coordinates is performed here.
     *  - Shear components in S use stress-strain convention (no engineering shear factor).
     *
     * Assumptions and validity:
     *  - Material symmetry is orthotropic in the principal directions.
     *  - Reciprocal Poisson relationships are enforced: n_yx = n_xy * Ey/Ex, etc., ensuring S is symmetric.
     *  - All moduli must be positive and physically consistent (E > 0, G > 0, |nu| < 0.5 for isotropy,
     *    and thermodynamic admissibility for orthotropy).
     *
     * Output:
     *  - S: On return, contains the 6x6 compliance matrix S in principal coordinates.
     */
    virtual void ComputeCompliance3DPrincipal(TPZFMatrix<REAL> &S) const;

    /**
     * @brief Compute compliance matrix for 3D elasticity
     * @param S Output compliance matrix
     * this method will rotate the compliance matrix from principal to global coordinates
     */
    virtual void ComputeComplianceMatrix3D(TPZFMatrix<REAL> &S) const;
    /**
     * @brief Compute compliance matrix S = C^{-1}
     * @param S Output compliance matrix
     *
     * If dimension==2, behavior depends on fPlaneStress flag:
     *  - fPlaneStress == true  -> plane stress
     *  - fPlaneStress == false -> plane strain
     */
    void ComputeComplianceMatrix(TPZFMatrix<REAL> &S) const;

    /**
     * @brief Set plane stress/strain behavior for 2D problems
     * @param planeStress true -> plane stress, false -> plane strain (plane elasticity)
     */
    void SetPlaneStress(bool planeStress) { fPlaneStress = planeStress; }

    /**
     * @brief Query plane stress flag (only meaningful when dimension==2)
     * @return true if plane stress, false if plane strain
     */
    bool IsPlaneStress() const { return fPlaneStress; }
    
    /**
     * @brief Compute strain-displacement matrix B
     * @param gradPhi Gradient of shape functions
     * @param B Output strain-displacement matrix
     */
    virtual void ComputeStrainDisplacementMatrix(const TPZFMatrix<REAL> &gradPhi, 
                                                TPZFMatrix<REAL> &B) const;
    
    /// Check if material is isotropic
    bool IsIsotropic() const { return fIsIsotropic; }
    
    /// Get problem dimension (2 or 3)
    int Dimension() const { return fDimension; }

    /// Set problem dimension; only 2 or 3 are allowed
    void SetDimension(int dim) {
#ifdef PZDEBUG
        if (dim != 2 && dim != 3) {
            DebugStop();
            return;
        }
#endif
        fDimension = dim;
    }

    /**
     * @brief Compute stress from strain using constitutive law
     * @param strain Strain vector in Voigt notation
     * @param stress Output stress vector in Voigt notation
     */
    virtual void ComputeStress(const TPZFMatrix<REAL> &strain, TPZFMatrix<REAL> &stress) const;

    /// @brief Compute div (sigma) as a function of the second derivatives of the displacement field
    /// @param d2udx2 Second derivatives of displacement field
    /// @param divsigma Output div(sigma) vector
    virtual void ComputeDivSigma(const TPZVec<REAL> &x, const TPZVec<TPZFNMatrix<9,REAL> > &d2udx2, TPZVec<REAL> &divsigma) const;

    /// @brief Compute sigma_z component from 2D strain
    /// @param strain in 2 dimensions
    /// @return computed sigma_z for plane strain/ returns 0 for plane stress
    /// throws an exception if the dimension is not 2
    /// this method assumes all components of young modulus and poisson ratio are set
    REAL ComputeSigmaZ(TPZFMatrix<REAL> &strain2D) const;

    /// @brief Compute epsilon_z component from 2D strain
    /// @param strain in 2 dimensions
    /// throws an exception if the dimension is not 2
    /// this method compute the deformation in z for plane stress/ returns 0 for plane strain
    REAL ComputeStrainZ(TPZFMatrix<REAL> &strain2D) const;
    /**
     * @brief Set principal axis directions for orthotropic material
     * @param axis1 First principal direction
     * @param axis2 Second principal direction
     * @param axis3 Third principal direction
     */
    void SetPrincipalAxes(const TPZVec<REAL> &axis1, 
                         const TPZVec<REAL> &axis2, 
                         const TPZVec<REAL> &axis3);
    
    /// Set rotation matrix directly (overrides principal axes)
    /// equivalent to setting fPrincipalAxes where the axes are stored column wise
    void SetRotationMatrix(const TPZFMatrix<REAL> &R);
    
    /// Compute rotation matrix from principal to global coordinates
    void ComputeRotationMatrix(TPZFMatrix<REAL> &R, bool transpose) const;

    /// @brief Set a 2D rotation angle (in radians) for principal axes
    /// @param angle Rotation angle in radians
    void SetPrincipalAxes2D(REAL angle);

    /// Build rotation matrix in Voigt notation
    // the matrix T is such that: {e_global} = T {e_principal}
    // the method assumes R is a proper three dimensional rotation matrix
    // T is 6x6
    // if S is a stress tensor associated with system B 
    // R transforms a vector in system A to a vector in system B
    // then T(R) S(B) = R S R^T is a stress tensor associated with system A
    static void BuildRotationVoigt3D(const TPZFMatrix<REAL> &R, TPZFMatrix<REAL> &T, bool isshearstrain = false);

    /// Build rotation matrix in Voigt notation
    // the matrix T is such that: {e_global} = T {e_principal}
    // the method assumes R is a proper three dimensional rotation matrix which is assumed to represent
    // a rotation in the plane
    // T is 3x3
    // if S is a stress tensor associated with system B 
    // R transforms a vector in system A to a vector in system B
    // then T(R) S(B) = R S R^T is a stress tensor associated with system A
    static void BuildRotationVoigt2D(const TPZFMatrix<REAL> &R, TPZFMatrix<REAL> &T, bool isshearstrain = false);

    /// @brief Rotate stress vector from principal to global coordinates
    /// @param tensor_principal : tensor in principal directions
    /// @param tensor_global : tensor in global directions
    /// @param transpose : if true, transpose the rotation (in this case, rotate from global to principal)
    /// the rotation is the one defined by fPrincipalAxes (RP) (principal to global)
    /// if(transpose = false) tensor_global = T(RP) tensor_principal
    /// if(transpose = true) tensor_principal = T(RP)^T tensor_global
    void RotateVoigtVector(const TPZFMatrix<REAL> &tensor_principal, TPZFMatrix<REAL> &tensor_global, bool transpose, bool isshearstrain = false) const;
protected:
    
    int fDimension = 2; // Dimension of the problem (2D or 3D)

    /// Flag indicating if material is isotropic
    bool fIsIsotropic;
    
    /// Orthotropic elastic moduli
    REAL fEx, fEy, fEz;
    
    /// Orthotropic Poisson's ratios
    REAL fNuxy, fNuyz, fNuxz;
    
    /// Orthotropic shear moduli
    REAL fGxy, fGyz, fGxz;
    
    /// Principal axis directions (stored as column vectors)
    /// this matrix represents a rotation that takes vectors from principal to global coordinates
    TPZFNMatrix<9, REAL> fPrincipalAxes;


    /// Plane stress flag (only for 2D problems)
    bool fPlaneStress;
    
};

#endif // TPZLINEARELASTICITYCONSTITUTIVE_H