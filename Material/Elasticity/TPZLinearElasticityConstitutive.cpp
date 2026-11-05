#include "TPZLinearElasticityConstitutive.h"
#include "pzmatrix.h"
#include "pzvec.h"
#include <cmath>
#include <stdexcept>

TPZLinearElasticityConstitutive::TPZLinearElasticityConstitutive()
: fIsIsotropic(true),
    fEx(0.), fEy(0.), fEz(0.),
    fNuxy(0.), fNuyz(0.), fNuxz(0.),
    fGxy(0.), fGyz(0.), fGxz(0.),
    fPrincipalAxes(3,3,0.)
{
        // Default principal axes = identity
        fPrincipalAxes.Identity();
}

TPZLinearElasticityConstitutive::~TPZLinearElasticityConstitutive() = default;

TPZLinearElasticityConstitutive::TPZLinearElasticityConstitutive(REAL E, REAL nu)
: TPZLinearElasticityConstitutive()
{
        SetIsotropicProperties(E, nu);
}

TPZLinearElasticityConstitutive::TPZLinearElasticityConstitutive(REAL Ex, REAL Ey, REAL Ez,
                                                                                                                                 REAL nuxy, REAL nuyz, REAL nuxz,
                                                                                                                                 REAL Gxy, REAL Gyz, REAL Gxz)
: TPZLinearElasticityConstitutive()
{
        SetOrthotropicProperties(Ex, Ey, Ez, nuxy, nuyz, nuxz, Gxy, Gyz, Gxz);
}

    /// @brief Constructor for orthotropic material in 2D
    /// @param Ex Young's modulus in x direction
    /// @param Ey Young's modulus in y direction
    /// @param nuxy Poisson's ratio xy
    /// @param Gxy Shear modulus xy
    TPZLinearElasticityConstitutive::TPZLinearElasticityConstitutive(REAL Ex, REAL Ey,
                                   REAL nuxy,
                                   REAL Gxy, bool planeStress) 
                                   : fDimension(2), fIsIsotropic(false),
                                   fEx(Ex), fEy(Ey), fEz(0.),
                                   fNuxy(nuxy), fNuyz(0.), fNuxz(0.),
                                   fGxy(Gxy), fGyz(0.), fGxz(0.), fPlaneStress(planeStress), fPrincipalAxes(3,3,0.)
{
        // Default principal axes = identity
        fPrincipalAxes.Identity();
                                   }

TPZLinearElasticityConstitutive::TPZLinearElasticityConstitutive(const TPZLinearElasticityConstitutive &copy)
: fIsIsotropic(copy.fIsIsotropic),
    fEx(copy.fEx), fEy(copy.fEy), fEz(copy.fEz),
    fNuxy(copy.fNuxy), fNuyz(copy.fNuyz), fNuxz(copy.fNuxz),
    fGxy(copy.fGxy), fGyz(copy.fGyz), fGxz(copy.fGxz),
    fPrincipalAxes(copy.fPrincipalAxes)
{
}

TPZLinearElasticityConstitutive &
TPZLinearElasticityConstitutive::operator=(const TPZLinearElasticityConstitutive &copy)
{
        if (this != &copy) {
                fIsIsotropic = copy.fIsIsotropic;
                fEx = copy.fEx; fEy = copy.fEy; fEz = copy.fEz;
                fNuxy = copy.fNuxy; fNuyz = copy.fNuyz; fNuxz = copy.fNuxz;
                fGxy = copy.fGxy; fGyz = copy.fGyz; fGxz = copy.fGxz;
                fPrincipalAxes = copy.fPrincipalAxes;
                fPlaneStress = copy.fPlaneStress;
        }
        return *this;
}

void TPZLinearElasticityConstitutive::SetIsotropicProperties(REAL E, REAL nu)
{
        if (E <= 0.) {
                throw std::invalid_argument("Isotropic E must be positive");
        }
        if (nu <= -1.0 || nu >= 0.5) {
                throw std::invalid_argument("Isotropic nu must be in (-1, 0.5)");
        }
        fIsIsotropic = true;

        fEx = E; fEy = E; fEz = E;
        fNuxy = nu; fNuyz = nu; fNuxz = nu;

        const REAL G = E / (2. * (1. + nu));
        fGxy = G; fGyz = G; fGxz = G;

        // principal axes = identity by default
        fPrincipalAxes.Redim(3,3);
        fPrincipalAxes.Zero();
        fPrincipalAxes(0,0) = 1.; fPrincipalAxes(1,1) = 1.; fPrincipalAxes(2,2) = 1.;
}

void TPZLinearElasticityConstitutive::SetOrthotropicProperties(REAL Ex, REAL Ey, REAL Ez,
                                        REAL nuxy, REAL nuyz, REAL nuxz,
                                        REAL Gxy, REAL Gyz, REAL Gxz)
{
        if (Ex <= 0. || Ey <= 0) {
                throw std::invalid_argument("Orthotropic E's must be positive");
        }
        if (fDimension == 3 && Ez <= 0.) {
                throw std::invalid_argument("Orthotropic Ez must be positive in 3D");
        }   
        if (Gxy <= 0) {
                throw std::invalid_argument("Orthotropic shear moduli must be positive");
        }
        if (fDimension == 3 && (Gyz <= 0. || Gxz <= 0.)) {
                throw std::invalid_argument("Orthotropic Gyz and Gxz must be positive in 3D");
        }
        // Basic sanity for Poisson ratios
        if (std::abs(nuxy) >= 0.99 || std::abs(nuyz) >= 0.99 || std::abs(nuxz) >= 0.99) {
                throw std::invalid_argument("Orthotropic Poisson ratios must be in (-1,1)");
        }
        if(fDimension != 2 && fDimension !=3) {
            fDimension = 3; // Default to 3D if not set
        }
        fIsIsotropic = false;
        fEx = Ex; fEy = Ey; fEz = Ez;
        fNuxy = nuxy; fNuyz = nuyz; fNuxz = nuxz;
        fGxy = Gxy; fGyz = Gyz; fGxz = Gxz;

        // principal axes keep current value; user can set later
        if (fPrincipalAxes.Rows() != 3 || fPrincipalAxes.Cols() != 3) {
                fPrincipalAxes.Redim(3,3);
                fPrincipalAxes.Zero();
                fPrincipalAxes(0,0) = 1.; fPrincipalAxes(1,1) = 1.; fPrincipalAxes(2,2) = 1.;
        }
}

/**
 * @brief Set orthotropic material properties
 */
void TPZLinearElasticityConstitutive::SetOrthotropicProperties(REAL Ex, REAL Ey,
                                REAL nuxy,
                                REAL Gxy) {

    if (Ex <= 0. || Ey <= 0) {
            throw std::invalid_argument("Orthotropic E's must be positive");
    }
    if (Gxy <= 0) {
            throw std::invalid_argument("Orthotropic shear moduli must be positive");
    }
    // Basic sanity for Poisson ratios
    if (std::abs(nuxy) >= 0.99) {
            throw std::invalid_argument("Orthotropic Poisson ratios must be in (-1,1)");
    }
    fIsIsotropic = false;
    fDimension = 2;
    fPlaneStress = true; // Default to plane stress for 2D orthotropic
    fEx = Ex; fEy = Ey; fEz = Ex; // Plane stress/strain: Ez not used
    fNuxy = nuxy; fNuyz = nuxy; fNuxz = nuxy;
    fGxy = Gxy; fGyz = Gxy; fGxz = Gxy;
    // principal axes keep current value; user can set later
    if (fPrincipalAxes.Rows() != 3 || fPrincipalAxes.Cols() != 3) {
            fPrincipalAxes.Redim(3,3);
            fPrincipalAxes.Identity();
    }
}

void TPZLinearElasticityConstitutive::SetPrincipalAxes(const TPZVec<REAL> &axis1,
     const TPZVec<REAL> &axis2,
    const TPZVec<REAL> &axis3) {
        if (axis1.size() < 3 || axis2.size() < 3 || axis3.size() < 3) {
                throw std::invalid_argument("Principal axes must be 3D vectors");
        }
        // Normalize and assign as columns
        auto norm = [](REAL x, REAL y, REAL z) {
                return std::sqrt(x*x + y*y + z*z);
        };
        REAL n1 = norm(axis1[0], axis1[1], axis1[2]);
        REAL n2 = norm(axis2[0], axis2[1], axis2[2]);
        REAL n3 = norm(axis3[0], axis3[1], axis3[2]);
        if (n1 == 0. || n2 == 0. || n3 == 0.) {
                throw std::invalid_argument("Principal axes must be non-zero");
        }
        // Check orthogonality of the normalized principal axes
        const REAL tol = 1e-8;
        const REAL dot12 = (axis1[0]*axis2[0] + axis1[1]*axis2[1] + axis1[2]*axis2[2])/(n1*n2);
        const REAL dot13 = (axis1[0]*axis3[0] + axis1[1]*axis3[1] + axis1[2]*axis3[2])/(n1*n3);
        const REAL dot23 = (axis2[0]*axis3[0] + axis2[1]*axis3[1] + axis2[2]*axis3[2])/(n2*n3);
        if (std::abs(dot12) > tol || std::abs(dot13) > tol || std::abs(dot23) > tol) {
            throw std::invalid_argument("Principal axes must be mutually orthogonal (dot products not near zero)");
        }

        fPrincipalAxes.Redim(3,3);
        fPrincipalAxes(0,0) = axis1[0]/n1; fPrincipalAxes(1,0) = axis1[1]/n1; fPrincipalAxes(2,0) = axis1[2]/n1;
        fPrincipalAxes(0,1) = axis2[0]/n2; fPrincipalAxes(1,1) = axis2[1]/n2; fPrincipalAxes(2,1) = axis2[2]/n2;
        fPrincipalAxes(0,2) = axis3[0]/n3; fPrincipalAxes(1,2) = axis3[1]/n3; fPrincipalAxes(2,2) = axis3[2]/n3;
}

    /// @brief Set a 2D rotation angle (in radians) for principal axes
    /// @param angle Rotation angle in radians
    void TPZLinearElasticityConstitutive::SetPrincipalAxes2D(REAL angle)
    {
        // Compute the rotation matrix for the given angle
        TPZFMatrix<REAL> R(2,2);
        R(0,0) = std::cos(angle); R(0,1) = -std::sin(angle);
        R(1,0) = std::sin(angle); R(1,1) = std::cos(angle);

        // Set the rotation matrix
        SetRotationMatrix(R);
    }

void TPZLinearElasticityConstitutive::SetRotationMatrix(const TPZFMatrix<REAL> &R)
{
    // Accept 2D rotation by embedding it in 3D with unit z-axis
    if (R.Rows() == 2 && R.Cols() == 2) {
        TPZFNMatrix<9,REAL> R3(3,3,0.);
        R3(0,0) = R(0,0); R3(0,1) = R(0,1); R3(0,2) = 0.;
        R3(1,0) = R(1,0); R3(1,1) = R(1,1); R3(1,2) = 0.;
        R3(2,0) = 0.;     R3(2,1) = 0.;     R3(2,2) = 1.;
        fPrincipalAxes = R3;
        return;
    } else if (R.Rows() != 3 || R.Cols() != 3) {
            throw std::invalid_argument("Rotation matrix must be 3x3");
    }
    fPrincipalAxes = R;
}

    /// @brief Rotate stress vector from principal to global coordinates
    /// @param stress_original 
    /// @param stress_rotated 
    void TPZLinearElasticityConstitutive::RotateVoigtVector(const TPZFMatrix<REAL> &stress_original, TPZFMatrix<REAL> &stress_rotated, bool transpose, bool is_shear_strain) const
{
        // Rotate the stress vector using the rotation matrix
        TPZFMatrix<REAL> T;
        if(fDimension == 2) T.Redim(3,3);
        else T.Redim(6,6);

        if(transpose == false) {
            if(fDimension == 2) {
                BuildRotationVoigt2D(fPrincipalAxes, T, is_shear_strain);
            } else {
                BuildRotationVoigt3D(fPrincipalAxes, T, is_shear_strain);
            }
        }

        else {
            if(fDimension == 2) {
                TPZFNMatrix<3,REAL> R_transpose(3,3);
                ComputeRotationMatrix(R_transpose, true);
                BuildRotationVoigt2D(R_transpose, T, is_shear_strain);
            } else {
                TPZFNMatrix<9,REAL> R_transpose(3,3);
                ComputeRotationMatrix(R_transpose, true);
                BuildRotationVoigt3D(R_transpose, T, is_shear_strain);
            }
        }
        stress_rotated = T * stress_original;
}

void TPZLinearElasticityConstitutive::ComputeRotationMatrix(TPZFMatrix<REAL> &R, bool transpose) const
{
        // R maps principal (material) to global coordinates
        // Columns of fPrincipalAxes are unit vectors in global coords
        R.Redim(3,3);
        for (int i=0; i<3; i++) {
                for (int j=0; j<3; j++) {
                    if(!transpose) {
                        R(i,j) = fPrincipalAxes(i,j);
                    } else {
                        R(i,j) = fPrincipalAxes(j,i);
                    }
                }
        }
        // Optionally, could orthonormalize, but assume user provided orthonormal basis.
}

static void BuildIsotropicC(REAL E, REAL nu, TPZFMatrix<REAL> &C)
{
        const REAL lambda = (nu*E)/((1.+nu)*(1.-2.*nu));
        const REAL mu = E/(2.*(1.+nu));

        C.Redim(6,6);
        C.Zero();

        // Normal components
        C(0,0) = lambda + 2.*mu;
        C(1,1) = lambda + 2.*mu;
        C(2,2) = lambda + 2.*mu;

        C(0,1) = lambda; C(0,2) = lambda;
        C(1,0) = lambda; C(1,2) = lambda;
        C(2,0) = lambda; C(2,1) = lambda;

        // Shear components (engineering shear)
        C(3,3) = mu;
        C(4,4) = mu;
        C(5,5) = mu;
}
void TPZLinearElasticityConstitutive::ComputeCompliance3DPrincipal(TPZFMatrix<REAL> &S) const
{
        // Compliance in principal coordinates (engineering strains)
        S.Zero();

        // Reciprocal Poisson ratios for symmetric S:
        // nu_yx = nu_xy * Ey/Ex, etc.

        const REAL n_yx = fNuxy * fEy / fEx;
        const REAL n_zy = fNuyz * fEz / fEy;
        const REAL n_zx = fNuxz * fEz / fEx;

        S(0,0) = 1./fEx;
        S(1,1) = 1./fEy;
        S(2,2) = 1./fEz;

        S(0,1) = -n_yx/fEy;
        S(1,0) = -fNuxy/fEx;

        S(1,2) = -n_zy/fEz;
        S(2,1) = -fNuyz/fEy;

        S(0,2) = -n_zx/fEz;
        S(2,0) = -fNuxz/fEx;

        S(3,3) = 1./fGyz;
        S(4,4) = 1./fGxz;
        S(5,5) = 1./fGxy;
}

void TPZLinearElasticityConstitutive::ComputeStiffnessMatrix3D(TPZFMatrix<REAL> &C) const
{
        if (fIsIsotropic) {
                BuildIsotropicC(fEx, fNuxy, C);
                return;
        }

        // Build C in principal coordinates
        TPZFNMatrix<36,REAL> C_princ;
        ComputeStiffnessMatrix3DPrincipal(C_princ);
        // Rotate to global: C_global = T * C_princ * T^T
        TPZFNMatrix<9,REAL> R(3,3), RT(3,3);
        ComputeRotationMatrix(R, false);
        ComputeRotationMatrix(RT, true);
        TPZFNMatrix<36,REAL> T, TT;
        BuildRotationVoigt3D(R, T, false);
        BuildRotationVoigt3D(RT, TT, true);

        TPZFNMatrix<36,REAL> tmp(6,6,0.), Cglob(6,6,0.);
        T.Multiply(C_princ, tmp, 0);          // tmp = T * C_princ
        tmp.Multiply(TT, Cglob, 0);            // Cglob = tmp * T^T  (A.Mult(B,result,1) means A*B^T)
        C = Cglob;
}


void TPZLinearElasticityConstitutive::ComputeComplianceMatrix3D(TPZFMatrix<REAL> &S) const
{
    // First compute compliance in principal coordinates
    TPZFMatrix<REAL> S_principal(6,6);
    ComputeCompliance3DPrincipal(S_principal);

    // Rotate compliance to global coordinates: S_global = T * S_principal * T^T
    TPZFMatrix<REAL> T(6,6), TT(6,6);
    TPZFNMatrix<9,REAL> R(3,3), R_transpose(3,3);
    ComputeRotationMatrix(R, false);
    ComputeRotationMatrix(R_transpose, true);
    BuildRotationVoigt3D(R, T, true); // true = shear strain
    BuildRotationVoigt3D(R_transpose, TT, false); // false = stress

    TPZFMatrix<REAL> temp(6,6);
    T.Multiply(S_principal, temp);
    temp.Multiply(TT, S);
}

#include <array>
#include <cmath>


void TPZLinearElasticityConstitutive::BuildRotationVoigt2D(const TPZFMatrix<REAL> &Rot, TPZFMatrix<REAL> &T, bool isshearstrain)
{
    double c = Rot(0,0);
    double s = Rot(1,0);

    // Rotation matrix in Voigt notation
    double R[3][3] = {
        { c*c, s*s, 2*c*s },
        { s*s, c*c, -2*c*s },
        { -c*s, c*s, c*c - s*s }
    };

    T.Redim(3,3);
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            T(i,j) = R[i][j];
        }
    }
    if(isshearstrain) {
        // Adjust for engineering shear strain factors
        for(int i=0; i<2; i++) {
            for(int j=2; j<3; j++) {
                T(i,j) *= 0.5;
                T(j,i) *= 2.0;
            }
        }
    }

}
void TPZLinearElasticityConstitutive::BuildRotationVoigt3D(const TPZFMatrix<REAL> &R, TPZFMatrix<REAL> &T, bool isshearstrain)
{
        // Build 6x6 transformation matrix for engineering strain Voigt vector:
        // e = [exx, eyy, ezz, 2*eyz, 2*exz, 2*exy]^T
        // e_global = T * e_principal
        T.Redim(6,6);
        T.Zero();

        const REAL r11=R(0,0), r12=R(0,1), r13=R(0,2);
        const REAL r21=R(1,0), r22=R(1,1), r23=R(1,2);
        const REAL r31=R(2,0), r32=R(2,1), r33=R(2,2);

        // Helper to fill T using standard tensor transformation rules
        auto Q = [&](int i,int j){ return R(i,j); };

        // Normal components
        T(0,0) = Q(0,0)*Q(0,0); T(0,1) = Q(0,1)*Q(0,1); T(0,2) = Q(0,2)*Q(0,2);
        T(0,3) = Q(0,1)*Q(0,2); T(0,4) = Q(0,0)*Q(0,2); T(0,5) = Q(0,0)*Q(0,1);

        T(1,0) = Q(1,0)*Q(1,0); T(1,1) = Q(1,1)*Q(1,1); T(1,2) = Q(1,2)*Q(1,2);
        T(1,3) = Q(1,1)*Q(1,2); T(1,4) = Q(1,0)*Q(1,2); T(1,5) = Q(1,0)*Q(1,1);

        T(2,0) = Q(2,0)*Q(2,0); T(2,1) = Q(2,1)*Q(2,1); T(2,2) = Q(2,2)*Q(2,2);
        T(2,3) = Q(2,1)*Q(2,2); T(2,4) = Q(2,0)*Q(2,2); T(2,5) = Q(2,0)*Q(2,1);

        // Shear components (engineering shear: 2*e_ij)
        // Row for 2*e_yz
        T(3,0) = 2.*Q(1,0)*Q(2,0);
        T(3,1) = 2.*Q(1,1)*Q(2,1);
        T(3,2) = 2.*Q(1,2)*Q(2,2);
        T(3,3) = Q(1,1)*Q(2,2)+Q(1,2)*Q(2,1);
        T(3,4) = Q(1,0)*Q(2,2)+Q(1,2)*Q(2,0);
        T(3,5) = Q(1,0)*Q(2,1)+Q(1,1)*Q(2,0);

        // Row for 2*e_xz
        T(4,0) = 2.*Q(0,0)*Q(2,0);
        T(4,1) = 2.*Q(0,1)*Q(2,1);
        T(4,2) = 2.*Q(0,2)*Q(2,2);
        T(4,3) = Q(0,1)*Q(2,2)+Q(0,2)*Q(2,1);
        T(4,4) = Q(0,0)*Q(2,2)+Q(0,2)*Q(2,0);
        T(4,5) = Q(0,0)*Q(2,1)+Q(0,1)*Q(2,0);

        // Row for 2*e_xy
        T(5,0) = 2.*Q(0,0)*Q(1,0);
        T(5,1) = 2.*Q(0,1)*Q(1,1);
        T(5,2) = 2.*Q(0,2)*Q(1,2);
        T(5,3) = Q(0,1)*Q(1,2)+Q(0,2)*Q(1,1);
        T(5,4) = Q(0,0)*Q(1,2)+Q(0,2)*Q(1,0);
        T(5,5) = Q(0,0)*Q(1,1)+Q(0,1)*Q(1,0);
    T.Transpose();
    if(isshearstrain) {
        // Adjust for engineering shear strain factors
        for(int i=0; i<3; i++) {
            for(int j=3; j<6; j++) {
                T(i,j) *= 0.5;
                T(j,i) *= 2.0;
            }
        }
    }
}

void TPZLinearElasticityConstitutive::ComputeStiffnessMatrix(TPZFMatrix<REAL> &C) const
{
    if (fDimension == 2) {
        if (C.Rows() != 3 || C.Cols() != 3) {
            DebugStop();
        }
        C.Zero();
    } else if (fDimension == 3) {
        if (C.Rows() != 6 || C.Cols() != 6) {
            DebugStop();
        }
        C.Zero();
    } else {
        throw std::invalid_argument("Unsupported dimension for constitutive matrix.");
    }
    if(fDimension == 2) {
        if(fPlaneStress) {
            ComputeStiffnessMatrixPlaneStress(C);
        }
        else {
            ComputeStiffnessMatrixPlaneStrain(C);
        }
    }
    else if(fDimension == 3) {
        ComputeStiffnessMatrix3D(C);
    }
    else {
        throw std::invalid_argument("Unsupported dimension for stress computation.");
    }
}

void TPZLinearElasticityConstitutive::ComputeStrainDisplacementMatrix(const TPZFMatrix<REAL> &gradPhi,
                                                                                                                                            TPZFMatrix<REAL> &B) const
{
        // gradPhi: 3 x nshape, gradients in global coordinates
        const int nrow = gradPhi.Rows();
        const int nshape = gradPhi.Cols();
        if (nrow != 3) {
                throw std::invalid_argument("gradPhi must be 3 x nshape");
        }
        B.Redim(6, 3*nshape);
        B.Zero();

        for (int a = 0; a < nshape; a++) {
                const REAL dphidx = gradPhi(0,a);
                const REAL dphidy = gradPhi(1,a);
                const REAL dphidz = gradPhi(2,a);

                const int col = 3*a;

                // Normal strains
                B(0, col+0) = dphidx; // exx
                B(1, col+1) = dphidy; // eyy
                B(2, col+2) = dphidz; // ezz

                // Engineering shear strains
                B(5, col+0) = dphidy; // 2*exy -> dphi/dy in u_x
                B(5, col+1) = dphidx; // and dphi/dx in u_y

                B(4, col+0) = dphidz; // 2*exz -> dphi/dz in u_x
                B(4, col+2) = dphidx; // and dphi/dx in u_z

                B(3, col+1) = dphidz; // 2*eyz -> dphi/dz in u_y
                B(3, col+2) = dphidy; // and dphi/dy in u_z
        }
}

// When asked for your name
void TPZLinearElasticityConstitutive::ComputeStiffnessMatrixPlaneStress(TPZFMatrix<REAL> &C) const
{
    if (fIsIsotropic) {
        const REAL mu = fEx / (2. * (1. + fNuxy));

        C.Redim(3, 3);
        C.Zero();
        REAL denom = 1. - fNuxy * fNuxy;
        // Normal components
        C(0, 0) = fEx / denom;
        C(1, 1) = fEx / denom;
        C(0, 1) = fEx * fNuxy / denom;
        C(1, 0) = fEx * fNuxy / denom;

        // Shear component
        C(2, 2) = mu;
    } else {
    
        // Extract in-plane 3x3 for [exx, eyy, gamma_xy] -> indices [0,1,5]
        REAL nu_yx = fNuxy * fEy / fEx;
        REAL denom = 1. - fNuxy * nu_yx;
        C.Redim(3, 3);
        C(0, 0) = fEx / denom;
        C(0, 1) = nu_yx * fEx / denom;
        C(0, 2) = 0.0;

        C(1, 0) = fEy * fNuxy / denom;
        C(1, 1) = fEy / denom;
        C(1, 2) = 0.0;
        C(2, 2) = fGxy;
        // C.Print("C2D before rotation: ",std::cout);
        TPZFNMatrix<9,REAL> T(3,3), TT(3,3), R(2,2), RT(2,2),Ctemp(3,3);
        for(int i=0; i<2; i++) {
            for(int j=0; j<2; j++) {
                RT(i,j) = fPrincipalAxes(i,j);
                R(i,j) = fPrincipalAxes(j,i);
            }
        }
        BuildRotationVoigt2D(R,T,false);
        BuildRotationVoigt2D(RT,TT,true);
        T.Multiply(C,Ctemp);
        Ctemp.Multiply(TT,C);

    }
}

void TPZLinearElasticityConstitutive::ComputeStiffnessMatrixPlaneStrain(TPZFMatrix<REAL> &C) const
{
    if (fIsIsotropic) {
        const REAL E = fEx;
        const REAL nu = fNuxy;
        const REAL factor = E / ((1. + nu) * (1. - 2. * nu));

        C.Redim(3, 3);
        C.Zero();

        // Normal components
        C(0, 0) = factor * (1. - nu);
        C(1, 1) = factor * (1. - nu);
        C(0, 1) = factor * nu;
        C(1, 0) = factor * nu;

        // Shear component
        C(2, 2) = factor * (1. - 2. * nu) / 2.;
    } else {
        // For orthotropic materials, compute full 6x6 matrix then extract plane strain components
        TPZFNMatrix<36,REAL> C_full;
        ComputeStiffnessMatrix3D(C_full);

//        C_full.Print("C3D full: ",std::cout);
        // For plane strain: ezz = 0, gamma_xz = 0, gamma_yz = 0
        // Extract the 3x3 block corresponding to [exx, eyy, gamma_xy]
        C.Redim(3, 3);
        C(0, 0) = C_full(0, 0); // exx-exx
        C(0, 1) = C_full(0, 1); // exx-eyy
        C(1, 0) = C_full(1, 0); // eyy-exx
        C(1, 1) = C_full(1, 1); // eyy-eyy
        C(0, 2) = C_full(0, 5); // exx-gamma_xy
        C(1, 2) = C_full(1, 5); // eyy-gamma_xy
        C(2, 0) = C_full(5, 0); // gamma_xy-exx
        C(2, 1) = C_full(5, 1); // gamma_xy-eyy
        C(2, 2) = C_full(5, 5); // gamma_xy-gamma_xy
    }
}

void TPZLinearElasticityConstitutive::ComputeStress(const TPZFMatrix<REAL> &strain, TPZFMatrix<REAL> &stress) const {
    int dimension = fDimension; // Use fDimension to determine the problem's dimension

    if (strain.Rows() != (dimension == 3 ? 6 : 3)) {
        throw std::invalid_argument("Strain vector must have 6 components for 3D or 3 components for 2D.");
    }

    if (stress.Rows() != (dimension == 3 ? 6 : 3)) {
        stress.Redim(dimension == 3 ? 6 : 3, 1);
    }

    TPZFNMatrix<36,REAL> C;
    if (dimension == 2) {
        if (fPlaneStress) {
            ComputeStiffnessMatrixPlaneStress(C);
        } else {
            ComputeStiffnessMatrixPlaneStrain(C);
        }
    } else if (dimension == 3) {
        ComputeStiffnessMatrix3D(C);
    } else {
        throw std::invalid_argument("Unsupported dimension for stress computation.");
    }
    // Stress = C * Strain
    C.Multiply(strain, stress);
}

    /**
     * @brief Compute compliance matrix S = C^{-1}
     * @param S Output compliance matrix
     *
     * If dimension==2, behavior depends on fPlaneStress flag:
     *  - fPlaneStress == true  -> plane stress
     *  - fPlaneStress == false -> plane strain
     */
void TPZLinearElasticityConstitutive::ComputeComplianceMatrix(TPZFMatrix<REAL> &S) const {
    int dimension = fDimension;
    if (dimension == 2) {
        if (S.Rows() != 3 || S.Cols() != 3) {
            S.Redim(3,3);
        }
    } else if (dimension == 3) {
        if (S.Rows() != 6 || S.Cols() != 6) {
            S.Redim(6,6);
        }
    } else {
        throw std::invalid_argument("Unsupported dimension for compliance computation.");
    }

    if (dimension == 2) {
        if (fPlaneStress) {
            TPZFNMatrix<36,REAL> S3D(6,6);
            ComputeComplianceMatrix3D(S3D);
            // Extract 2D plane stress compliance from 3D compliance
            S(0,0) = S3D(0,0); S(0,1) = S3D(0,1); S(0,2) = S3D(0,5);
            S(1,0) = S3D(1,0); S(1,1) = S3D(1,1); S(1,2) = S3D(1,5);
            S(2,0) = S3D(5,0); S(2,1) = S3D(5,1); S(2,2) = S3D(5,5);
        } else if (!fPlaneStress) {
            DebugStop(); // Plane stress compliance not implemented
        } else {
            DebugStop(); // Plane strain compliance not implemented
        }
    } else if (dimension == 3) {
        ComputeComplianceMatrix3D(S);
    } else {
        throw std::invalid_argument("Unsupported dimension for compliance computation.");
    }
}

void TPZLinearElasticityConstitutive::ComputeStiffnessMatrix3DPrincipal(TPZFMatrix<REAL> &C) const {

    C.Redim(6, 6);
    C.Zero();

    // Reciprocity conditions
    double nu21 = fNuxy * fEy / fEx;
    double nu32 = fNuyz * fEz / fEy;
    double nu31 = fNuxz * fEz / fEx;

    // Denominator Δ
    double Delta =
        1.0 - fNuxy*nu21 - fNuyz*nu32 - nu31*fNuxz
        - 2.0*fNuxy*fNuyz*nu31;

    // --- Normal terms ---

    C(0,0) =  (1 - fNuyz*nu32) * fEx / Delta;
    C(0,1) =  (nu21 + fNuyz*nu31) * fEx / Delta;
    C(0,2) =  (nu31 + nu21*nu32) * fEx / Delta;

    C(1,0) =  (fNuxy + fNuxz*nu32) * fEy / Delta;
    C(1,1) =  (1 - fNuxz*nu31)    * fEy / Delta;
    C(1,2) =  (nu32 + fNuxy*nu31) * fEy / Delta;

    C(2,0) =  (fNuxz + fNuxy*fNuyz) * fEz / Delta;
    C(2,1) =  (fNuyz + nu21*fNuxz) * fEz / Delta;
    C(2,2) =  (1 - fNuxy*nu21)    * fEz / Delta;

    // --- Zero normal–shear couplings ---
    for(int i = 0; i < 3; ++i)
        for(int j = 3; j < 6; ++j)
            C(i,j) = C(j,i) = 0.0;

    // --- Shear terms ---
    C(3,3) = fGyz;
    C(4,4) = fGxz;
    C(5,5) = fGxy;

    // Other shear off-diagonals are zero
    C(3,4) = C(3,5) = 0.0;
    C(4,3) = C(4,5) = 0.0;
    C(5,3) = C(5,4) = 0.0;

}

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
void TPZLinearElasticityConstitutive::TensorToVoigt(const TPZFMatrix<REAL> &T, TPZFMatrix<REAL> &v, bool engineeringShear) const {
    if (fDimension == 3) {
        if (T.Rows() != 3 || T.Cols() != 3) {
            throw std::invalid_argument("Input tensor T must be 3x3 for 3D.");
        }
        v.Redim(6, 1);
        v(0,0) = T(0,0); // xx
        v(1,0) = T(1,1); // yy
        v(2,0) = T(2,2); // zz
        v(3,0) = engineeringShear ? T(1,2)+T(2,1) : 0.5*(T(1,2)+T(2,1)); // yz
        v(4,0) = engineeringShear ? T(0,2)+T(2,0) : 0.5*(T(0,2)+T(2,0)); // xz
        v(5,0) = engineeringShear ? T(0,1)+T(1,0) : 0.5*(T(0,1)+T(1,0)); // xy
    } else if (fDimension == 2) {
        if (T.Rows() != 2 || T.Cols() != 2) {
            throw std::invalid_argument("Input tensor T must be 2x2 for 2D.");
        }
        v.Redim(3, 1);
        v(0,0) = T(0,0); // xx
        v(1,0) = T(1,1); // yy
        v(2,0) = engineeringShear ? T(0,1)+T(1,0) : 0.5*(T(0,1)+T(1,0)); // xy
    } else {
        throw std::invalid_argument("Unsupported dimension for TensorToVoigt.");
    }
}

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
void TPZLinearElasticityConstitutive::VoigtToTensor(const TPZFMatrix<REAL> &v, TPZFMatrix<REAL> &T, bool engineeringShear) const {
    if (fDimension == 3) {
        if (v.Rows() != 6 || v.Cols() != 1) {
            throw std::invalid_argument("Input Voigt vector v must be 6x1 for 3D.");
        }
        T.Redim(3, 3);
        T(0,0) = v(0,0); // xx
        T(1,1) = v(1,0); // yy
        T(2,2) = v(2,0); // zz
        T(1,2) = engineeringShear ? v(3,0)/2. : v(3,0); // yz
        T(2,1) = T(1,2);
        T(0,2) = engineeringShear ? v(4,0)/2. : v(4,0); // xz
        T(2,0) = T(0,2);
        T(0,1) = engineeringShear ? v(5,0)/2. : v(5,0); // xy
        T(1,0) = T(0,1);
    } else if (fDimension == 2) {
        if (v.Rows() != 3 || v.Cols() != 1) {
            throw std::invalid_argument("Input Voigt vector v must be 3x1 for 2D.");
        }
        T.Redim(2, 2);
        T(0,0) = v(0,0); // xx
        T(1,1) = v(1,0); // yy
        T(0,1) = engineeringShear ? v(2,0)/2. : v(2,0); // xy
        T(1,0) = T(0,1);
    } else {
        throw std::invalid_argument("Unsupported dimension for VoigtToTensor.");
    }
}

/// @brief Compute sigma_z component from 2D strain
/// @param strain in 2 dimensions
/// @return computed sigma_z for plane strain/ returns 0 for plane stress
/// throws an exception if the dimension is not 2
/// this method assumes all components of young modulus and poisson ratio are set
REAL TPZLinearElasticityConstitutive::ComputeSigmaZ(TPZFMatrix<REAL> &strain2D) const {
    if(fDimension != 2) {
        throw std::invalid_argument("ComputeSigmaZ only valid for 2D problems.");
    }
    if(strain2D.Rows() !=3 || strain2D.Cols() !=1) {
        throw std::invalid_argument("Input strain2D must be a 3x1 Voigt vector [exx, eyy, gamma_xy]^T.");
    }

    if(fPlaneStress) {
        // For plane stress, sigma_z = 0
        return 0.0;
    } else {
        // For plane strain, compute sigma_z
        TPZFNMatrix<36,REAL> C_full;
        ComputeStiffnessMatrix3D(C_full);

        REAL exx = strain2D(0,0);
        REAL eyy = strain2D(1,0);
        REAL ezz = 0.0;
        REAL gamma_xy = strain2D(2,0);
        REAL sigma_z = C_full(2,0)*exx + C_full(2,1)*eyy + C_full(2,2)*ezz + C_full(2,5)*gamma_xy;
        return sigma_z;
    }
}

/// @brief Compute epsilon_z component from 2D strain
/// @param strain in 2 dimensions
/// throws an exception if the dimension is not 2
/// this method compute the deformation in z for plane stress/ returns 0 for plane strain
REAL TPZLinearElasticityConstitutive::ComputeStrainZ(TPZFMatrix<REAL> &strain2D) const {
    if(fDimension != 2) {
        throw std::invalid_argument("ComputeStrainZ only valid for 2D problems.");
    }
    if(strain2D.Rows() !=3 || strain2D.Cols() !=1) {
        throw std::invalid_argument("Input strain2D must be a 3x1 Voigt vector [exx, eyy, gamma_xy]^T.");
    }

    if(fPlaneStress) {
        // For plane stress, compute epsilon_z
        if(fIsIsotropic == false) {
            DebugStop(); // Not implemented for anisotropic materials
        }
        REAL exx = strain2D(0,0);
        REAL eyy = strain2D(1,0);
//        std::cout << - fNuxy/(1-fNuxy) << " "<< exx << "+"<< - fNuxy/(1-fNuxy) << " "<< eyy <<std::endl;
        REAL ezz = - fNuxy/(1-fNuxy)*(exx + eyy);
        return ezz;
    } else {
        // For plane strain, epsilon_z = 0
        return 0.0;
    }
}   

/// @brief Compute div (sigma) as a function of the second derivatives of the displacement field
/// @param d2udx2 Second derivatives of displacement field
/// @param divsigma Output div(sigma) vector
void TPZLinearElasticityConstitutive::ComputeDivSigma(const TPZVec<REAL> &x, const TPZVec<TPZFNMatrix<9,REAL> > &d2udx2, TPZVec<REAL> &divsigma) const {
    int nstrains = (fDimension == 3) ? 6 : 3;
    TPZFNMatrix<36,REAL> D(nstrains,nstrains);
    ComputeStiffnessMatrix(D);
    TPZFNMatrix<18,REAL> gradstrain(nstrains,fDimension), gradsigma;
    for(int j=0; j<fDimension; j++) {
        gradstrain(0,j) = d2udx2[0](0,j);
        gradstrain(1,j) = d2udx2[1](1,j);
        if(fDimension == 2) {
            gradstrain(2,j) = d2udx2[0](1,j)+d2udx2[1](0,j);
        }
        if(fDimension ==3) {
            gradstrain(2,j) = d2udx2[2](2,j);
            gradstrain(3,j) = d2udx2[1](2,j)+d2udx2[2](1,j);
            gradstrain(4,j) = d2udx2[0](2,j)+d2udx2[2](0,j);
            gradstrain(5,j) = d2udx2[0](1,j)+d2udx2[1](0,j);
        }
    }
    D.Multiply(gradstrain,gradsigma);
    divsigma.resize(fDimension);
    if(fDimension == 2) {
        divsigma[0]=gradsigma(0,0) + gradsigma(2,1);
        divsigma[1]=gradsigma(1,1) + gradsigma(2,0);
    }
    if(fDimension == 3) {
        divsigma[0]=gradsigma(0,0) + gradsigma(4,2) + gradsigma(5,1);
        divsigma[1]=gradsigma(1,1) + gradsigma(3,2) + gradsigma(5,0);
        divsigma[2]=gradsigma(2,2) + gradsigma(3,1) + gradsigma(4,0);
    }
}
