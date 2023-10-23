/**
 * @file
 * @brief Contains implementations of the TPZMixedElasticityND methods.
 */

#include "TPZMixedElasticityND.h"
#include "pzelmat.h"
#include "TPZBndCondT.h"
#include "pzaxestools.h"
#include "TPZMatWithMem.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzerror.h"
#include "TPZMaterialDataT.h"
#include <math.h>

#include "pzlog.h"
//#include "meshgen.h"
#include "TPZAnalyticSolution.h"
#include "TPZLapackEigenSolver.h"
#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.mixedelasticity"));
#endif

#include <fstream>
using namespace std;

TPZMixedElasticityND::TPZMixedElasticityND() : TBase() {
    fE_const = -1.; // Young modulus
    fnu_const = -1.; // Poisson coefficient
    fForce[0] = 0.; // X component of the body force
    fForce[1] = 0.; // Y component of the body force
    fForce[2] = 0.; // Z component of the body force - not used for this class
    flambda_const = 0.;
    fmu_const = 0.;

    fPlaneStress = 0;
    fMatrixA = 0.;
}

TPZMixedElasticityND::TPZMixedElasticityND(int id, int dimension) : TBase(id), fDimension(dimension) {
    fE_const = -1.; // Young modulus
    fnu_const = -1.; // Poisson coefficient
    fForce[0] = 0.; // X component of the body force
    fForce[1] = 0.; // Y component of the body force
    fForce[2] = 0.; // Z component of the body force - not used for this class
    flambda_const = 0.;
    fmu_const = 0.;

    fPlaneStress = 0;
    fMatrixA = 0.;
}

TPZMixedElasticityND::TPZMixedElasticityND(int id, REAL E, REAL nu, REAL fx, REAL fy, int planestress, int dimension) : TBase(id), fDimension(dimension) {
    this->SetElasticity(E, nu);
    fForce[0] = fx; // X component of the body force
    fForce[1] = fy; // Y component of the body force
    fForce[2] = 0.;
    fPlaneStress = planestress;
}

TPZMixedElasticityND::~TPZMixedElasticityND() {
}

void TPZMixedElasticityND::FillDataRequirements(TPZVec<TPZMaterialDataT<STATE> > &datavec) const {
    int nref = datavec.size();
    for (int i = 0; i < nref; i++) {
        datavec[i].SetAllRequirements(false);
        datavec[i].fNeedsNeighborSol = false;
        datavec[i].fNeedsNeighborCenter = false;
        datavec[i].fNeedsNormal = false;
    }
}

int TPZMixedElasticityND::NStateVariables() const {
    return fDimension;
}

////////////////////////////////////////////////////////////////////

// Divergence on deformed element

void TPZMixedElasticityND::ComputeDivergenceOnDeformed(TPZVec<TPZMaterialDataT<STATE>> &datavec, TPZFMatrix<STATE> &DivergenceofPhi) {
    //itapopo conferir esse método. Foi copiado do TPZDarcyFlow3D
    int sigmaBlock = 0;

    // Getting test and basis functions
    //TPZFMatrix<REAL> phiSigmaH1 = datavec[sigmaBlock].phi; // For H1  test functions Q
    //TPZFMatrix<REAL> dphiSigmaH1 = datavec[sigmaBlock].dphi; // Derivative For H1  test functions
    TPZFMatrix<REAL> dphiSigmaH1axes = datavec[sigmaBlock].dphix; // Derivative For H1  test functions

    TPZFNMatrix<660> GradphiuH1;
    TPZAxesTools<REAL>::Axes2XYZ(dphiSigmaH1axes, GradphiuH1, datavec[sigmaBlock].axes);

    int nphiuHdiv = datavec[sigmaBlock].fVecShapeIndex.NElements();

    DivergenceofPhi.Resize(nphiuHdiv, 1);

    REAL JacobianDet = datavec[sigmaBlock].detjac;

    TPZFMatrix<REAL> Qaxes = datavec[sigmaBlock].axes;
    TPZFMatrix<REAL> QaxesT;
    TPZFMatrix<REAL> Jacobian = datavec[sigmaBlock].jacobian;
    TPZFMatrix<REAL> JacobianInverse = datavec[sigmaBlock].jacinv;

    TPZFMatrix<REAL> GradOfX;
    TPZFMatrix<REAL> GradOfXInverse;
    TPZFMatrix<REAL> VectorOnMaster;
    TPZFMatrix<REAL> VectorOnXYZ(3, 1, 0.0);
    Qaxes.Transpose(&QaxesT);
    QaxesT.Multiply(Jacobian, GradOfX);
    JacobianInverse.Multiply(Qaxes, GradOfXInverse);

    int ivectorindex = 0;
    int ishapeindex = 0;

    {
        for (int iq = 0; iq < nphiuHdiv; iq++) {
            ivectorindex = datavec[sigmaBlock].fVecShapeIndex[iq].first;
            ishapeindex = datavec[sigmaBlock].fVecShapeIndex[iq].second;

            VectorOnXYZ(0, 0) = datavec[sigmaBlock].fDeformedDirections(0, ivectorindex);
            VectorOnXYZ(1, 0) = datavec[sigmaBlock].fDeformedDirections(1, ivectorindex);
            VectorOnXYZ(2, 0) = datavec[sigmaBlock].fDeformedDirections(2, ivectorindex);

            GradOfXInverse.Multiply(VectorOnXYZ, VectorOnMaster);
            VectorOnMaster *= JacobianDet;

            /* Contravariant Piola mapping preserves the divergence */
            DivergenceofPhi(iq, 0) = (1.0 / JacobianDet) * (GradphiuH1(0, ishapeindex) * VectorOnMaster(0, 0) +
                    GradphiuH1(1, ishapeindex) * VectorOnMaster(1, 0) +
                    GradphiuH1(2, ishapeindex) * VectorOnMaster(2, 0));
        }
    }
    return;

}

void TPZMixedElasticityND::ElasticityModulusTensor(TPZFMatrix<STATE> &MatrixElast, TElasticityAtPoint &elast) {
    //Matrix modulus Voigt notation:
    int matdim = (fDimension == 3) ? 9 : 4;
    MatrixElast.Redim(matdim, matdim);
    REAL mu = elast.fmu;
    REAL lambda = elast.flambda;
    if(fDimension == 2)
    {
        if (!fPlaneStress) {
            // plane strain
            MatrixElast(Exx, Exx) = 1. / (2. * mu) - lambda / (2. * mu * (2. * lambda + 2. * mu));
            MatrixElast(Exx, Eyy) = -lambda / (2. * mu * (2. * lambda + 2. * mu));
            MatrixElast(Eyy, Exx) = MatrixElast(Exx, Eyy);
            MatrixElast(Eyy, Eyy) = MatrixElast(Exx, Exx);
            MatrixElast(Exy, Exy) = 1. / (2. * mu);
            MatrixElast(Eyx, Eyx) = MatrixElast(Exy, Exy);
        } else {
            // plane stress
//            MatrixElast(Exx, Exx) = (1 - elast.fnu*elast.fnu) / elast.fE;
//            MatrixElast(Exx, Eyy) = -elast.fnu * (1 + elast.fnu) / elast.fE;
//            MatrixElast(Eyy, Exx) = MatrixElast(Exx, Eyy);
//            MatrixElast(Eyy, Eyy) = MatrixElast(Exx, Exx);
//            MatrixElast(Exy, Exy) = 1. / (2. * mu);
//            MatrixElast(Eyx, Eyx) = MatrixElast(Exy, Exy);
            MatrixElast(Exx, Exx) = 1. / elast.fE;
            MatrixElast(Exx, Eyy) = -elast.fnu / elast.fE;
            MatrixElast(Eyy, Exx) = MatrixElast(Exx, Eyy);
            MatrixElast(Eyy, Eyy) = MatrixElast(Exx, Exx);
            MatrixElast(Exy, Exy) = (1.+elast.fnu) / elast.fE;
            MatrixElast(Eyx, Eyx) = MatrixElast(Exy, Exy);
        }
    } else if (fDimension == 3)
    {
        // sig = lambda tr(E) I + 2 mu E
        // E = 1/(2 mu) (sig - lambda /(3 lambda + 2 mu) tr(sig) I)
		// NOTE: NS double-checked and this is correct
        MatrixElast(Exx, Exx) = 1./(2.*mu)*(1.-lambda/(3.*lambda+2.*mu));
        MatrixElast(Exx, Eyy) = 1./(2.*mu)*(-lambda/(3.*lambda+2.*mu));
        MatrixElast(Exx, Ezz) = MatrixElast(Exx, Eyy);
        MatrixElast(Eyy, Eyy) = MatrixElast(Exx, Exx);
        MatrixElast(Eyy, Exx) = MatrixElast(Exx, Eyy);
        MatrixElast(Eyy, Ezz) = MatrixElast(Exx, Eyy);
        MatrixElast(Ezz, Ezz) = MatrixElast(Exx, Exx);
        MatrixElast(Ezz, Exx) = MatrixElast(Exx, Eyy);
        MatrixElast(Ezz, Eyy) = MatrixElast(Exx, Eyy);
        MatrixElast(Exy, Exy) = 1. / (2. * mu);
        MatrixElast(Exz, Exz) = MatrixElast(Exy, Exy);
        MatrixElast(Eyx, Eyx) = MatrixElast(Exy, Exy);
        MatrixElast(Eyz, Eyz) = MatrixElast(Exy, Exy);
        MatrixElast(Ezx, Ezx) = MatrixElast(Exy, Exy);
        MatrixElast(Ezy, Ezy) = MatrixElast(Exy, Exy);
    }
    else
    {
        DebugStop();
    }
#ifdef LOG4CXX
    if(logdata->isDebugEnabled()){
        std::stringstream sout;
        MatrixElast.Print("A",sout,EMathematicaInput);
        LOGPZ_DEBUG(logdata,sout.str());
    }
#endif //LOG4CXX
}

void TPZMixedElasticityND::ComputeDeformationVector(TPZVec<STATE> &PhiStress, TPZVec<STATE> &APhiStress, TElasticityAtPoint &elast) {
    int matdim = (fDimension == 3) ? 9 : 4;
    TPZFNMatrix<81, STATE> MatrixElast(matdim, matdim, 0.);
    ElasticityModulusTensor(MatrixElast, elast);
    for (int iq = 0; iq < matdim; iq++) {
        APhiStress[iq] = 0.;
        for (int jq = 0; jq < matdim; jq++) {
            APhiStress[iq] += MatrixElast(iq, jq) * PhiStress[jq];
        }
    }
}

void TPZMixedElasticityND::ComputeStressVector(TPZVec<STATE> &Deformation, TPZVec<STATE> &Stress, TElasticityAtPoint &elast) {
    if(fDimension == 2)
    {
        if (fPlaneStress) {
            Stress[Exx] = elast.fE / (1 - elast.fnu * elast.fnu)*(Deformation[Exx] + elast.fnu * Deformation[Eyy]);
            Stress[Eyy] = elast.fE / (1. - elast.fnu * elast.fnu)*(elast.fnu * Deformation[Exx] + Deformation[Eyy]);
            Stress[Exy] = elast.fE / (1. + elast.fnu) * Deformation[Exy];
            Stress[Eyx] = elast.fE / (1. + elast.fnu) * Deformation[Eyx];
        } else {
            Stress[Exx] = elast.fE / ((1 + elast.fnu)*(1 - 2. * elast.fnu))*((1 - elast.fnu) * Deformation[Exx] + elast.fnu * Deformation[Eyy]);
            Stress[Eyy] = elast.fE / ((1 + elast.fnu)*(1 - 2. * elast.fnu))*((1 - elast.fnu) * Deformation[Eyy] + elast.fnu * Deformation[Exx]);
            Stress[Exy] = elast.fE / (1. + elast.fnu) * Deformation[Exy];
            Stress[Eyx] = elast.fE / (1. + elast.fnu) * Deformation[Eyx];
        }
    } else if(fDimension == 3)
    {
        REAL mu = elast.fmu;
        REAL lambda = elast.flambda;
        STATE trE = Deformation[Exx]+Deformation[Eyy]+Deformation[Ezz];
        Stress[Exx] = lambda*trE+2.*mu*Deformation[Exx];
        Stress[Eyy] = lambda*trE+2.*mu*Deformation[Eyy];
        Stress[Ezz] = lambda*trE+2.*mu*Deformation[Ezz];
        Stress[Exy] = 2.*mu*Deformation[Exy];
        Stress[Exz] = 2.*mu*Deformation[Exz];
        Stress[Eyx] = 2.*mu*Deformation[Eyx];
        Stress[Eyz] = 2.*mu*Deformation[Eyz];
        Stress[Ezx] = 2.*mu*Deformation[Ezx];
        Stress[Ezy] = 2.*mu*Deformation[Ezy];
    }
}

void TPZMixedElasticityND::Print(std::ostream &out) const {
    out << "name of material : " << Name() << "\n";
    out << "properties : \n";
    out << "\tE_const   = " << fE_const << endl;
    out << "\tnu_const   = " << fnu_const << endl;
    out << "\tF   = " << fForce[0] << ' ' << fForce[1] << endl;
}

/**
 * Computes the spectral radius of the compliance tensor using Mandel notation.
 * Only the isotropic, minor-symmetric case is implemented.
 * A sigma = 1/(2mu) * (sigma - lambda/(n*lambda+2mu) * tr sigma I)
 * 
 * A[ijkl] = 1/(2mu) * (.5*(delta[ik]*delta[jl]+delta[il]*delta[kj]) - lambda/(n*lambda+2mu) * delta[ij]*delta[kl])
 * @param x
 * @return 
 */
CSTATE TPZMixedElasticityND::GetMaxComplianceEigenvalue(TPZVec<REAL> &x) const {
    if (fElasticity){
        // Not implemented for variable elasticity coefficients yet!
        
    // routines for the non-symmetric case
//    auto symsym = [&complianceElement](unsigned int i, unsigned int j, unsigned int k, unsigned int l) {
//        auto Aijkl = complianceElement(i,j,k,l);
//        auto Ajikl = complianceElement(j,i,k,l);
//        auto Aijlk = complianceElement(i,j,l,k);
//        auto Ajilk = complianceElement(j,i,l,k);
//        return 0.25*(Aijkl + Ajikl + Aijlk + Ajilk);
//    };
//    
//    auto symskw = [&complianceElement](unsigned int i, unsigned int j, unsigned int k, unsigned int l) {
//        auto Aijkl = complianceElement(i,j,k,l);
//        auto Ajikl = complianceElement(j,i,k,l);
//        auto Aijlk = complianceElement(i,j,l,k);
//        auto Ajilk = complianceElement(j,i,l,k);
//        return 0.25*(Aijkl + Ajikl - Aijlk - Ajilk);
//    };
//    
//    auto skwsym = [&complianceElement](unsigned int i, unsigned int j, unsigned int k, unsigned int l) {
//        auto Aijkl = complianceElement(i,j,k,l);
//        auto Ajikl = complianceElement(j,i,k,l);
//        auto Aijlk = complianceElement(i,j,l,k);
//        auto Ajilk = complianceElement(j,i,l,k);
//        return 0.25*(Aijkl - Ajikl + Aijlk - Ajilk);
//    };
//    
//    auto skwskw = [&complianceElement](unsigned int i, unsigned int j, unsigned int k, unsigned int l) {
//        auto Aijkl = complianceElement(i,j,k,l);
//        auto Ajikl = complianceElement(j,i,k,l);
//        auto Aijlk = complianceElement(i,j,l,k);
//        auto Ajilk = complianceElement(j,i,l,k);
//        return 0.25*(Aijkl - Ajikl - Aijlk + Ajilk);
//    };
        DebugStop();
    } else {
        TPZAutoPointer<TPZMatrix<STATE>> A(new TPZFNMatrix<36, STATE>(6,6,0.));

        auto delta = [](unsigned int i, unsigned int j) {
            return (i == j) ? 1. : 0.;
        };

        auto complianceElement = [&delta, this](unsigned int i, unsigned int j, unsigned int k, unsigned int l) {
            return (1. / (2 * fmu_const))*((delta(j, k) * delta(j, l) + delta(i, l) * delta(k, j)) / 2. - fmu_const / (fDimension * flambda_const + 2 * fmu_const) * delta(i, j) * delta(k, l));
        };

        // Assuming the 4th order tensor is minor-symmetric (Aijkl=Ajikl=Ajilk), which is the case for linear elasticity
        A->operator ()(0, 0) = complianceElement(1, 1, 1, 1);
        A->operator ()(0, 1) = complianceElement(1, 1, 2, 2);
        A->operator ()(0, 2) = complianceElement(1, 1, 3, 3);
        A->operator ()(0, 3) = complianceElement(1, 1, 2, 3) * M_SQRT2;
        A->operator ()(0, 4) = complianceElement(1, 1, 3, 1) * M_SQRT2;
        A->operator ()(0, 5) = complianceElement(1, 1, 1, 2) * M_SQRT2;
        A->operator ()(1, 0) = complianceElement(2, 2, 1, 1);
        A->operator ()(1, 1) = complianceElement(2, 2, 2, 2);
        A->operator ()(1, 2) = complianceElement(2, 2, 3, 3);
        A->operator ()(1, 3) = complianceElement(2, 2, 2, 3) * M_SQRT2;
        A->operator ()(1, 4) = complianceElement(2, 2, 3, 1) * M_SQRT2;
        A->operator ()(1, 5) = complianceElement(2, 2, 1, 2) * M_SQRT2;
        A->operator ()(2, 0) = complianceElement(3, 3, 1, 1);
        A->operator ()(2, 1) = complianceElement(3, 3, 2, 2);
        A->operator ()(2, 2) = complianceElement(3, 3, 3, 3);
        A->operator ()(2, 3) = complianceElement(3, 3, 2, 3) * M_SQRT2;
        A->operator ()(2, 4) = complianceElement(3, 3, 3, 1) * M_SQRT2;
        A->operator ()(2, 5) = complianceElement(3, 3, 1, 2) * M_SQRT2;
        A->operator ()(3, 0) = complianceElement(2, 3, 1, 1) * M_SQRT2;
        A->operator ()(3, 1) = complianceElement(2, 3, 2, 2) * M_SQRT2;
        A->operator ()(3, 2) = complianceElement(2, 3, 3, 3) * M_SQRT2;
        A->operator ()(3, 3) = complianceElement(2, 3, 2, 3)*2;
        A->operator ()(3, 4) = complianceElement(2, 3, 3, 1)*2;
        A->operator ()(3, 5) = complianceElement(2, 3, 1, 2)*2;
        A->operator ()(4, 0) = complianceElement(3, 1, 1, 1) * M_SQRT2;
        A->operator ()(4, 1) = complianceElement(3, 1, 2, 2) * M_SQRT2;
        A->operator ()(4, 2) = complianceElement(3, 1, 3, 3) * M_SQRT2;
        A->operator ()(4, 3) = complianceElement(3, 1, 2, 3)*2;
        A->operator ()(4, 4) = complianceElement(3, 1, 3, 1)*2;
        A->operator ()(4, 5) = complianceElement(3, 1, 1, 2)*2;
        A->operator ()(5, 0) = complianceElement(1, 2, 1, 1) * M_SQRT2;
        A->operator ()(5, 1) = complianceElement(1, 2, 2, 2) * M_SQRT2;
        A->operator ()(5, 2) = complianceElement(1, 2, 3, 3) * M_SQRT2;
        A->operator ()(5, 3) = complianceElement(1, 2, 2, 3)*2;
        A->operator ()(5, 4) = complianceElement(1, 2, 3, 1)*2;
        A->operator ()(5, 5) = complianceElement(1, 2, 1, 2)*2;
        
        TPZLapackEigenSolver<STATE> eigensolver;
        eigensolver.SetMatrixA(A);
        TPZManVector<CSTATE, 6> eigenvalues;
        eigensolver.Solve(eigenvalues);
        unsigned int maxnormindex = 0;
        STATE max_norm = norm(eigenvalues[maxnormindex]);
        for (unsigned int i = 1; i < 6; ++i) {
            STATE this_norm = norm(eigenvalues[i]);
            if (this_norm > max_norm){
                max_norm = this_norm;
                maxnormindex = i;
            }
        }
        return eigenvalues[maxnormindex];
    }
    return 0;
}

/// Transform a tensor to a Voigt notation

void TPZMixedElasticityND::ToVoigt(const TPZFMatrix<STATE> &S, TPZVec<STATE> &Svoigt) const {
    Svoigt[Exx] = S(0, 0);
    Svoigt[Exy] = S(0, 1);
    Svoigt[Eyx] = S(1, 0);
    Svoigt[Eyy] = S(1, 1);
    if (fDimension > 2) {
        Svoigt[Exz] = S(0, 2);
        Svoigt[Eyz] = S(1, 2);
        Svoigt[Ezx] = S(2, 0);
        Svoigt[Ezy] = S(2, 1);
        Svoigt[Ezz] = S(2, 2);
    }
}

/// Transform a Voigt notation to a tensor

void TPZMixedElasticityND::FromVoigt(const TPZVec<STATE> &Svoigt, TPZFMatrix<STATE> &S) const {
    S(0, 0) = Svoigt[Exx];
    S(0, 1) = Svoigt[Exy];
    S(1, 0) = Svoigt[Eyx];
    S(1, 1) = Svoigt[Eyy];
    if (fDimension > 2) {
        S(0, 2) = Svoigt[Exz];
        S(1, 2) = Svoigt[Eyz];
        S(2, 0) = Svoigt[Ezx];
        S(2, 1) = Svoigt[Ezy];
        S(2, 2) = Svoigt[Ezz];
    }
}

void TPZMixedElasticityND::MandelFromSecondOrderTensor(const TPZFMatrix<STATE> &S, TPZVec<STATE> &SMandel) const {
    SMandel[0] = S(0, 0); // S[1,1]
    SMandel[1] = S(1, 1); // S[2,2]
    if (fDimension == 2){
        SMandel[2] = M_SQRT2 * 0.5 * (S(0, 1)+S(1, 0)); // sqrt{2}/2 (S[1,2]+S[2,1])
        SMandel[3] = M_SQRT2 * 0.5 * (S(1, 0)-S(0, 1)); // sqrt{2}/2 (S[2,1]-S[1,2])
    } else {
        SMandel[2] = S(2, 2); // S[3,3]
        SMandel[3] = M_SQRT2 * 0.5 * (S(1, 2)+S(2, 1)); // sqrt{2}/2 (S[2,3]+S[3,2])
        SMandel[4] = M_SQRT2 * 0.5 * (S(2, 0)+S(0, 2)); // sqrt{2}/2 (S[3,1]+S[1,3])
        SMandel[5] = M_SQRT2 * 0.5 * (S(0, 1)+S(1, 0)); // sqrt{2}/2 (S[1,2]+S[2,1])
        SMandel[6] = M_SQRT2 * 0.5 * (S(2, 1)-S(1, 2)); // sqrt{2}/2 (S[3,2]-S[2,3])
        SMandel[7] = M_SQRT2 * 0.5 * (S(0, 2)-S(2, 0)); // sqrt{2}/2 (S[1,3]-S[3,1])
        SMandel[8] = M_SQRT2 * 0.5 * (S(1, 0)-S(0, 1)); // sqrt{2}/2 (S[2,1]-S[1,2])
    }
}

void TPZMixedElasticityND::SecondOrderTensorFromMandel(const TPZVec<STATE> &SMandel, TPZFMatrix<STATE> &S) const {
    S(0, 0) = SMandel[0];
    S(1, 1) = SMandel[1];
    if (fDimension == 2){
        S(0, 1) = M_SQRT1_2*(SMandel[2]-SMandel[3]);
        S(1, 0) = M_SQRT1_2*(SMandel[2]+SMandel[3]);
    } else {
        S(0, 1) = M_SQRT1_2*(SMandel[5]-SMandel[8]);
        S(1, 0) = M_SQRT1_2*(SMandel[5]+SMandel[8]);
        S(0, 2) = M_SQRT1_2*(SMandel[4]+SMandel[7]);
        S(1, 2) = M_SQRT1_2*(SMandel[3]-SMandel[6]);
        S(2, 0) = M_SQRT1_2*(SMandel[4]-SMandel[7]);
        S(2, 1) = M_SQRT1_2*(SMandel[3]+SMandel[6]);
        S(2, 2) = SMandel[2];
    }
}

void TPZMixedElasticityND::MandelFromSymmetricSecondOrderTensor(const TPZFMatrix<STATE> &S, TPZVec<STATE> &SMandel) const {
    SMandel[0] = S(0, 0); // S[1,1]
    SMandel[1] = S(1, 1); // S[2,2]
    if (fDimension == 2){
        SMandel[2] = M_SQRT2 * S(0, 1); // sqrt{2}/2 (S[1,2]+S[2,1])
    } else {
        SMandel[2] = S(2, 2); // S[3,3]
        SMandel[3] = M_SQRT2 * S(2, 1); // sqrt{2}/2 (S[2,3]+S[3,2])
        SMandel[4] = M_SQRT2 * S(2, 0); // sqrt{2}/2 (S[3,1]+S[1,3])
        SMandel[5] = M_SQRT2 * S(0, 1); // sqrt{2}/2 (S[1,2]+S[2,1])
    }
}

void TPZMixedElasticityND::SymmetricSecondOrderTensorFromMandel(const TPZVec<STATE> &SMandel, TPZFMatrix<STATE> &S) const {
    S(0, 0) = SMandel[0];
    S(1, 1) = SMandel[1];
    if (fDimension == 2){
        S(0, 1) = S(1, 0) = M_SQRT1_2 * SMandel[2];
    } else {
        S(2, 2) = SMandel[2];
        S(1, 2) = S(2, 1) = M_SQRT1_2 * SMandel[3];
        S(0, 2) = S(2, 0) = M_SQRT1_2 * SMandel[4];
        S(0, 1) = S(1, 0) = M_SQRT1_2 * SMandel[5];
    }
}

/** @brief Calculates the element stiffness matrix using 3 spaces - Stress tensor, displacement, and skew-symmetric tensor (for weak symmetry) */
void TPZMixedElasticityND::Contribute_3spaces(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    REAL R = datavec[0].x[0];
    // Setting the phi's
    // E
    TPZFMatrix<REAL> &phiS = datavec[0].phi;
    TPZFMatrix<REAL> &dphiS = datavec[0].dphix;
    // U
    TPZFMatrix<REAL> &phiU = datavec[1].phi;
//    TPZFMatrix<REAL> &dphiU = datavec[1].dphix;
    // P
    TPZFMatrix<REAL> &phiP = datavec[2].phi;

    TElasticityAtPoint elast(fE_const, fnu_const);
    if (fElasticity) {
        TPZManVector<STATE, 3> result(2);
        TPZFNMatrix<4, STATE> Dres(0, 0);
        fElasticity(datavec[0].x, result, Dres);
        REAL E = result[0];
        REAL nu = result[1];
        elast = TElasticityAtPoint(E, nu);
    }
    //    
//    TPZFNMatrix<220, REAL> dphiSx(fDimension, dphiS.Cols());
//    TPZAxesTools<REAL>::Axes2XYZ(dphiS, dphiSx, datavec[0].axes);

//    TPZFNMatrix<220, REAL> dphiUx(fDimension, phiU.Cols());
//    TPZAxesTools<REAL>::Axes2XYZ(dphiU, dphiUx, datavec[1].axes);

//    TPZFNMatrix<220, REAL> dphiPx(fDimension, phiP.Cols());
//    TPZAxesTools<REAL>::Axes2XYZ(dphiU, dphiPx, datavec[2].axes);

    int nshapeS, nshapeU, nshapeP;
    nshapeS = datavec[0].fVecShapeIndex.NElements();
    nshapeU = datavec[1].phi.Rows();
    nshapeP = datavec[2].phi.Rows();
    const int firstequation_S = 0;
    const int firstequation_U = firstequation_S + nshapeS*fDimension;
    //const int firstequation_P = firstequation_U + nshapeU*fDimension;
    
    // number of asymetric tensors for each shape function
    int nrotations = (fDimension == 3) ? 3 : 1;

    TPZManVector<STATE, 3>  force(fDimension, 0.);
    
    TPZFNMatrix<3,STATE>    phiSi(3, 1, 0.), 
                            phiSj(fDimension, 1, 0.);
    
    TPZManVector<STATE, 3>  divSi1x(fDimension, 0.), 
                            divSi1y(fDimension, 0.), 
                            divSi1z(fDimension, 0.);
    
    int voigtdim = fDimension*fDimension;
    TPZManVector<STATE, 9>  phiSi1x(voigtdim, 0.0), phiSi1y(voigtdim, 0.0), phiSi1z(voigtdim, 0.0),
                            phiSj1x(voigtdim, 0.0), phiSj1y(voigtdim, 0.0), phiSj1z(voigtdim, 0.0),
                            phiPj1x(voigtdim, 0.0), phiPj1y(voigtdim, 0.0), phiPj1z(voigtdim, 0.0);
    
    TPZManVector<STATE, 3>  phiUi(fDimension, 0.0), 
                            phiUj(fDimension, 0.0), 
                            phiUj1x(fDimension, 0.0), 
                            phiUj1y(fDimension, 0.0), 
                            phiUj1z(fDimension, 0.0);;
    
    TPZFNMatrix<3,STATE>    phiPi(fDimension, 1, 0.0), 
                            phiPj(fDimension, 1, 0.0);
    
    TPZFNMatrix<3, REAL>    ivecS(3, 1, 0.);

    force = fForce;
    if (this->HasForcingFunction()) {
        fForcingFunction(datavec[0].x, force);
#ifdef LOG4CXX
        if (logdata->isDebugEnabled()) {
            std::stringstream sout;
            sout << " x = " << datavec[0].x << " force = " << force << std::endl;
            LOGPZ_DEBUG(logdata, sout.str())
        }
#endif
    }
//    datavec[0].ComputeFunctionDivergence();
    for (int i = 0; i < nshapeS; i++) {
        int iphi = datavec[0].fVecShapeIndex[i].second;
        int ivec = datavec[0].fVecShapeIndex[i].first;

        REAL divSi = 0.;

        for (int e = 0; e < 3; e++) {
            ivecS(e, 0) = datavec[0].fDeformedDirections(e, ivec);
            phiSi(e, 0) = phiS(iphi, 0) * ivecS(e, 0);
        }
//        TPZFNMatrix<3, REAL> axesvec(fDimension, 1, 0.);
//        datavec[0].axes.Multiply(ivecS, axesvec);
//
//        //calculando div(Si)
//        for (int f = 0; f < fDimension; f++) {
//            divSi += axesvec(f, 0) * dphiS(f, iphi);
//        }
        divSi = datavec[0].divphi(i);

        TPZFNMatrix<9, STATE> phiTensx(fDimension, fDimension, 0.), phiTensy(fDimension, fDimension, 0.),
        phiTensz(fDimension,fDimension,0.);

        for (int d = 0; d < fDimension; d++) {
            phiTensx(0, d) = phiSi(d, 0);
            phiTensy(1, d) = phiSi(d, 0);
            if (fDimension == 3) {
                phiTensz(2, d) = phiSi(d, 0);
            }
        }
        ToVoigt(phiTensx, phiSi1x);
        ToVoigt(phiTensy, phiSi1y);
        if (fDimension == 3) {
            ToVoigt(phiTensz, phiSi1z);
        }

        divSi1x[0] = divSi;
        divSi1y[1] = divSi;
        if(fDimension == 3) divSi1z[2] = divSi;

        // matrix K11 - (Matrix A * stress tensor) x test-function stress tensor
        for (int j = 0; j < nshapeS; j++) {
            int jphi = datavec[0].fVecShapeIndex[j].second;
            int jvec = datavec[0].fVecShapeIndex[j].first;

            for (int e = 0; e < fDimension; e++) {
                phiSj(e, 0) = phiS(jphi, 0) * datavec[0].fDeformedDirections(e, jvec);
            }

            
            TPZFNMatrix<9, STATE> phjTensx(fDimension, fDimension, 0.), phjTensy(fDimension, fDimension, 0.)
            , phjTensz(fDimension, fDimension, 0.);

            for(int d = 0; d< fDimension; d++)
            {
                phjTensx(0,d) = phiSj(d,0);
                phjTensy(1,d) = phiSj(d,0);
                if(fDimension == 3)
                {
                    phjTensz(2,d) = phiSj(d,0);
                }
            }
            ToVoigt(phjTensx, phiSj1x);
            ToVoigt(phjTensy, phiSj1y);
            if(fDimension == 3)
            {
                ToVoigt(phjTensz, phiSj1z);
            }

            //Multiply by Lamé parameters
            TPZManVector<STATE, 9> AphiSi1x(voigtdim, 0.0), AphiSi1y(voigtdim, 0.0), AphiSi1z(voigtdim, 0.0);

            ComputeDeformationVector(phiSi1x, AphiSi1x, elast);
            ComputeDeformationVector(phiSi1y, AphiSi1y, elast);
            if(fDimension == 3)
            {
                ComputeDeformationVector(phiSi1z, AphiSi1z, elast);
            }

            STATE valxx = InnerVec(AphiSi1x, phiSj1x);
            STATE valxy = InnerVec(AphiSi1x, phiSj1y);
            STATE valyx = InnerVec(AphiSi1y, phiSj1x);
            STATE valyy = InnerVec(AphiSi1y, phiSj1y);

            //Matrix K11
            if (fAxisSymmetric) {
                valxx /= R;
                valxy /= R;
                valyx /= R;
                valyy /= R;
            }
            ek(fDimension * i, fDimension * j) += weight * valxx;
            ek(fDimension * i, fDimension * j + 1) += weight * valxy;
            ek(fDimension * i + 1, fDimension * j) += weight * valyx;
            ek(fDimension * i + 1, fDimension * j + 1) += weight * valyy;
            if(fDimension == 3)
            {
                STATE valxz = InnerVec(AphiSi1x, phiSj1z);
                ek(fDimension * i, fDimension * j+2) += weight * valxz;
                STATE valyz = InnerVec(AphiSi1y, phiSj1z);
                ek(fDimension * i+1, fDimension * j+2) += weight * valyz;
                STATE valzx = InnerVec(AphiSi1z, phiSj1x);
                ek(fDimension * i + 2, fDimension * j) += weight * valzx;
                STATE valzy = InnerVec(AphiSi1z, phiSj1y);
                ek(fDimension * i + 2, fDimension * j + 1) += weight * valzy;
                STATE valzz = InnerVec(AphiSi1z, phiSj1z);
                ek(fDimension * i + 2, fDimension * j + 2) += weight * valzz;

            }
        }

        // matrix K21 and K12 - divergent of test-function stress tensor * displacement vector
        for (int j = 0; j < nshapeU; j++) {
            phiUj1x[0] = phiU(j, 0);
            phiUj1y[1] = phiU(j, 0);
            
            STATE valx = weight * InnerVec(divSi1x, phiUj1x);
            STATE valy = weight * InnerVec(divSi1y, phiUj1y);

            //position K21
            ek(fDimension * j + nshapeS * fDimension, fDimension * i) += valx;
            ek(fDimension * j + 1 + nshapeS * fDimension, fDimension * i + 1) += valy;

            //position K12
            ek(fDimension * i, fDimension * j + nshapeS * fDimension) += valx;
            ek(fDimension * i + 1, fDimension * j + 1 + nshapeS * fDimension) += valy;
            
            if(fDimension == 3)
            {
                phiUj1z[2] = phiU(j,0);
                REAL valz = weight * InnerVec(divSi1z, phiUj1z);
                ek(fDimension * j + 2 + nshapeS * fDimension, fDimension * i + 2) += valz;
                ek(fDimension * i + 2, fDimension * j + 2 + nshapeS * fDimension) += valz;
            }

        }

        // matrix K31 and K13 - test-function stress tensor x rotation tensor p
        for (int j = 0; j < nshapeP; j++) {
            TPZFNMatrix<9, STATE> phiPTensor(fDimension, fDimension, 0.);
            phiPTensor(0, 1) = phiP(j, 0);
            phiPTensor(1, 0) = -phiP(j, 0);
            ToVoigt(phiPTensor, phiPj1x);

            STATE valxx = InnerVec(phiSi1x, phiPj1x);
            STATE valxy = InnerVec(phiSi1y, phiPj1x);
            //Matrix K31
            ek(nrotations*j + nshapeS * fDimension + nshapeU * fDimension, fDimension * i) += weight * valxx;
            ek(nrotations*j + nshapeS * fDimension + nshapeU * fDimension, fDimension * i + 1) += weight * valxy;

            //Matrix K13
            ek(fDimension * i, nrotations*j + nshapeS * fDimension + nshapeU * fDimension) += weight * valxx;
            ek(fDimension * i + 1, nrotations*j + nshapeS * fDimension + nshapeU * fDimension) += weight * valxy;
            
            if(fDimension == 3)
            {
                STATE valxz = InnerVec(phiSi1z, phiPj1x);
                ek(nrotations*j + nshapeS * fDimension + nshapeU * fDimension, fDimension * i + 2) += weight * valxz;
                ek(fDimension * i + 2, nrotations*j + nshapeS * fDimension + nshapeU * fDimension) += weight * valxz;
                for(int d=0; d<2; d++)
                {
                    phiPTensor.Zero();
                    phiPTensor(d, 2) = phiP(j, 0);
                    phiPTensor(2, d) = -phiP(j, 0);
                    ToVoigt(phiPTensor, phiPj1x);

                    STATE valxx = InnerVec(phiSi1x, phiPj1x);
                    STATE valxy = InnerVec(phiSi1y, phiPj1x);
                    STATE valxz = InnerVec(phiSi1z, phiPj1x);

                    //Matrix K31
                    ek(nrotations*j + d + 1 + nshapeS * fDimension + nshapeU * fDimension, fDimension * i) += weight * valxx;
                    ek(nrotations*j + d + 1 + nshapeS * fDimension + nshapeU * fDimension, fDimension * i + 1) += weight * valxy;
                    ek(nrotations*j + d + 1 + nshapeS * fDimension + nshapeU * fDimension, fDimension * i + 2) += weight * valxz;

                    //Matrix K13
                    ek(fDimension * i, nrotations*j + d + 1 + nshapeS * fDimension + nshapeU * fDimension) += weight * valxx;
                    ek(fDimension * i + 1, nrotations*j + d + 1 + nshapeS * fDimension + nshapeU * fDimension) += weight * valxy;
                    ek(fDimension * i + 2, nrotations*j + d + 1 + nshapeS * fDimension + nshapeU * fDimension) += weight * valxz;

                }
            }
        }
    }

    for (int i = 0; i < nshapeU; i++) {
        phiUj1x[0] = phiU(i, 0);
        phiUj1y[1] = phiU(i, 0);

        // Load vector f:
        STATE factfx = -weight * phiUj1x[0] * force[0];
        STATE factfy = -weight * phiUj1y[1] * force[1];
        if (fAxisSymmetric) {
            factfx *= R;
            factfy *= R;
        }
        //if(factfx != 0) DebugStop();
        ef(nshapeS * fDimension + fDimension * i, 0) += factfx;
        ef(nshapeS * fDimension + fDimension * i + 1, 0) += factfy;
        
        if(fDimension == 3)
        {
            phiUj1z[2] = phiU(i, 0);
            STATE factfz = -weight * phiUj1z[2] * force[2];
            ef(nshapeS * fDimension + fDimension * i + 2, 0) += factfz;

        }

        //if(ef(nshapeS*2+2*j) != 0.) DebugStop();
        if (fAxisSymmetric)
        {
            for (int j = 0; j < nshapeU; j++) {
                ek(nshapeS * 2 + 2 * i, nshapeS * 2 + 2 * j) -= weight * phiU(i, 0) * phiU(j, 0) / (elast.fE * R);
            }
        }
    }
}

/** @brief Calculates the element stiffness matrix using 5 spaces - 3 from Contribute_3spaces() + Rigid body motions, and distributed forces */
void TPZMixedElasticityND::Contribute_5spaces(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
    int nspaces = datavec.size();
    Contribute_3spaces(datavec, weight, ek, ef);

    if(nspaces == 3) return;
//    ofstream out1("matrix1.txt");
//    ek.Print("matrix1",out1,EFixedColumn);
//    out1.flush();

    // Distributed forces
    TPZFMatrix<REAL>& phiFRB = datavec[3].phi;
    const int nshapeFRB = phiFRB.Rows();
    // Rigid body displacements
    TPZFMatrix<REAL>& phiURB = datavec[4].phi;
    const int nshapeURB = phiURB.Rows();
    // U
    TPZFMatrix<REAL> &phiU = datavec[1].phi;
    

    const int nshapeS = datavec[0].fVecShapeIndex.NElements();
    const int nshapeU = datavec[1].phi.Rows();
    const int nshapeP = datavec[2].phi.Rows();
    const int nrotations = fDimension == 3 ? 3:1;
    const int firstequation_S = 0;
    const int firstequation_U = firstequation_S + nshapeS*fDimension;
    const int firstequation_P = firstequation_U + nshapeU*fDimension;
    int firstequation_FRB = firstequation_P + nshapeP*nrotations;
    
    ContributeRigidBodyMode(datavec, weight, ek, ef, firstequation_FRB, 3);
    if(nspaces == 5) return;
    if(fDimension == 2)
    {
        firstequation_FRB += 3;
    }
    else
    {
        firstequation_FRB += 12;
    }
    ContributeRigidBodyMode(datavec, weight, ek, ef, firstequation_FRB, 5);

}
    // contribute a distributed flux and average displacement/rotation
void TPZMixedElasticityND::ContributeRigidBodyMode(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, int first_eq, int RB_space)
{
    int nspaces = datavec.size();
    // Distributed forces
    TPZFMatrix<REAL>& phiFRB = datavec[RB_space].phi;
    const int nshapeFRB = phiFRB.Rows();
    // Rigid body displacements
    TPZFMatrix<REAL>& phiURB = datavec[RB_space+1].phi;
    const int nshapeURB = phiURB.Rows();
    // U
    TPZFMatrix<REAL> &phiU = datavec[1].phi;
    

    const int nshapeS = datavec[0].fVecShapeIndex.NElements();
    const int nshapeU = datavec[1].phi.Rows();
    const int nshapeP = datavec[2].phi.Rows();
    const int nrotations = fDimension == 3 ? 3:1;
    const int firstequation_S = 0;
    const int firstequation_U = firstequation_S + nshapeS*fDimension;
    const int firstequation_P = firstequation_U + nshapeU*fDimension;
    const int firstequation_FRB = first_eq;

    const int ncomponents = fDimension+nrotations;
    const int firstequation_URB = firstequation_FRB + nshapeFRB*ncomponents;

    if(RB_space == nspaces-2)
    {
        bool sizetest = ek.Cols() == firstequation_FRB + 2*ncomponents;
        if(!sizetest) DebugStop();
    }
    // get distance of integration point to center of element

    TPZManVector<REAL,3> xcenter = datavec[RB_space].XCenter;
    TPZManVector<REAL,3> x = datavec[0].x;
    TPZManVector<REAL,3> xparam = datavec[RB_space].xParametric;
    TPZManVector<REAL,3> delx(3,0.);
                            delx[0] = datavec[0].x[0] - xcenter[0];
                            delx[1] = datavec[0].x[1] - xcenter[1];
    if(fDimension==3){  delx[2] = datavec[0].x[2] - xcenter[2];}

    // Matrix K42 and K24 - Distributed forces x displacement
    for(int j=0; j<nshapeU; j++){
        // Matrix K24
        // Translation forces x displacement
        STATE val_f = weight*phiU(j,0);
        ek(0 + j*fDimension + firstequation_U, 0 + firstequation_FRB) += val_f;
        ek(1 + j*fDimension + firstequation_U, 1 + firstequation_FRB) += val_f;
        if(fDimension==3) 
        {
            ek(2 + j*fDimension + firstequation_U, 2 + firstequation_FRB) += val_f;
        }
        // Rotation forces  x displacement
        STATE val_mx = weight*phiU(j,0)*delx[0];
        STATE val_my = weight*phiU(j,0)*delx[1];
        STATE val_mz = weight*phiU(j,0)*delx[2];
        if(fDimension==3){
            ek(0 + j*fDimension + firstequation_U, 1 + fDimension + firstequation_FRB) -= val_mz;
            ek(1 + j*fDimension + firstequation_U, 0 + fDimension + firstequation_FRB) += val_mz;
            ek(1 + j*fDimension + firstequation_U, 2 + fDimension + firstequation_FRB) -= val_mx;
            ek(2 + j*fDimension + firstequation_U, 1 + fDimension + firstequation_FRB) += val_mx;
            ek(0 + j*fDimension + firstequation_U, 2 + fDimension + firstequation_FRB) += val_my;
            ek(2 + j*fDimension + firstequation_U, 0 + fDimension + firstequation_FRB) -= val_my;
        }else{
            ek(0 + j*fDimension + firstequation_U, 0 + fDimension + firstequation_FRB) -= val_my;
            ek(1 + j*fDimension + firstequation_U, 0 + fDimension + firstequation_FRB) += val_mx;
        }
        // Matrix K42
        // Translation forces
        ek(0 + firstequation_FRB, 0 + j*fDimension + firstequation_U) += val_f;
        ek(1 + firstequation_FRB, 1 + j*fDimension + firstequation_U) += val_f;
        if(fDimension==3) 
            {ek(2 + firstequation_FRB, 2 + j*fDimension + firstequation_U) += val_f;}
        // Rotation forces
        if(fDimension==3){
            ek(1 + fDimension + firstequation_FRB, 0 + j*fDimension + firstequation_U) -= val_mz;
            ek(0 + fDimension + firstequation_FRB, 1 + j*fDimension + firstequation_U) += val_mz;
            ek(2 + fDimension + firstequation_FRB, 1 + j*fDimension + firstequation_U) -= val_mx;
            ek(1 + fDimension + firstequation_FRB, 2 + j*fDimension + firstequation_U) += val_mx;
            ek(2 + fDimension + firstequation_FRB, 0 + j*fDimension + firstequation_U) += val_my;
            ek(0 + fDimension + firstequation_FRB, 2 + j*fDimension + firstequation_U) -= val_my;
        }else{
            ek(0 + fDimension + firstequation_FRB, 0 + j*fDimension + firstequation_U) -= val_my;
            ek(0 + fDimension + firstequation_FRB, 1 + j*fDimension + firstequation_U) += val_mx;
        }
    }
    // Matrix K45 and K54 - Distributed forces x rigid body displacement
    //  Matrix 45
    // Translations x Translations
    for(int d=0; d<fDimension; d++)
        {ek(d + firstequation_FRB, d + firstequation_URB) -= weight;}
    // Translations x Rotations
    STATE val_x = weight*delx[0];
    STATE val_y = weight*delx[1];
    STATE val_z = weight*delx[2];
    if(fDimension==3){
        ek(2 + firstequation_FRB, 1 + firstequation_URB + fDimension) += val_x;
        ek(1 + firstequation_FRB, 2 + firstequation_URB + fDimension) -= val_x;
        ek(0 + firstequation_FRB, 2 + firstequation_URB + fDimension) += val_y;
        ek(2 + firstequation_FRB, 0 + firstequation_URB + fDimension) -= val_y;
        ek(1 + firstequation_FRB, 0 + firstequation_URB + fDimension) += val_z;
        ek(0 + firstequation_FRB, 1 + firstequation_URB + fDimension) -= val_z;
    }else{
        ek(0 + firstequation_FRB, firstequation_URB + fDimension) += val_y;
        ek(1 + firstequation_FRB, firstequation_URB + fDimension) -= val_x;
    }
    // Rotations x Translations
    if(fDimension==3){
        ek(1 + firstequation_FRB + fDimension, 2 + firstequation_URB) += val_x;
        ek(2 + firstequation_FRB + fDimension, 1 + firstequation_URB) -= val_x;
        ek(2 + firstequation_FRB + fDimension, 0 + firstequation_URB) += val_y;
        ek(0 + firstequation_FRB + fDimension, 2 + firstequation_URB) -= val_y;
        ek(0 + firstequation_FRB + fDimension, 1 + firstequation_URB) += val_z;
        ek(1 + firstequation_FRB + fDimension, 0 + firstequation_URB) -= val_z;
    }else{
        ek(0 + firstequation_FRB + fDimension, 0 + firstequation_URB) += val_y;
        ek(0 + firstequation_FRB + fDimension, 1 + firstequation_URB) -= val_x;
    }
    // Rotations x Rotations
    STATE val2_xy = weight*(delx[0]*delx[0] + delx[1]*delx[1]);
    STATE val2_xz = weight*(delx[0]*delx[0] + delx[2]*delx[2]);
    STATE val2_yz = weight*(delx[1]*delx[1] + delx[2]*delx[2]);
    STATE val_xy = weight*(delx[0]*delx[1]);
    STATE val_xz = weight*(delx[0]*delx[2]);
    STATE val_yz = weight*(delx[1]*delx[2]);
    
    
    if(fDimension==3){
        ek(0 + firstequation_FRB + fDimension, 0 + firstequation_URB + fDimension) += val2_yz; 
        ek(0 + firstequation_FRB + fDimension, 1 + firstequation_URB + fDimension) -= val_xy; 
        ek(0 + firstequation_FRB + fDimension, 2 + firstequation_URB + fDimension) -= val_xz; 
        ek(1 + firstequation_FRB + fDimension, 0 + firstequation_URB + fDimension) -= val_xy; 
        ek(1 + firstequation_FRB + fDimension, 1 + firstequation_URB + fDimension) += val2_xz; 
        ek(1 + firstequation_FRB + fDimension, 2 + firstequation_URB + fDimension) -= val_yz; 
        ek(2 + firstequation_FRB + fDimension, 0 + firstequation_URB + fDimension) -= val_xz; 
        ek(2 + firstequation_FRB + fDimension, 1 + firstequation_URB + fDimension) -= val_yz; 
        ek(2 + firstequation_FRB + fDimension, 2 + firstequation_URB + fDimension) += val2_xy; 
    }else{
        ek(0 + firstequation_FRB + fDimension, 0 + firstequation_URB + fDimension) -= val2_xy; 
    }
    //  Matrix 54
    // Translations x Translations
    for(int d=0; d<fDimension; d++)
        {ek(d + firstequation_URB, d + firstequation_FRB) -= weight;}
    // Translations x Rotations
    if(fDimension==3){
        ek(1 + firstequation_URB + fDimension, 2 + firstequation_FRB) += val_x;
        ek(2 + firstequation_URB + fDimension, 1 + firstequation_FRB) -= val_x;
        ek(2 + firstequation_URB + fDimension, 0 + firstequation_FRB) += val_y;
        ek(0 + firstequation_URB + fDimension, 2 + firstequation_FRB) -= val_y;
        ek(0 + firstequation_URB + fDimension, 1 + firstequation_FRB) += val_z;
        ek(1 + firstequation_URB + fDimension, 0 + firstequation_FRB) -= val_z;
    }else{
        ek(firstequation_URB + fDimension, 0 + firstequation_FRB) += val_y;
        ek(firstequation_URB + fDimension, 1 + firstequation_FRB) -= val_x;
    }
    // Rotations x Translations
    if(fDimension==3){
        ek(2 + firstequation_URB, 1 + firstequation_FRB + fDimension) += val_x;
        ek(1 + firstequation_URB, 2 + firstequation_FRB + fDimension) -= val_x;
        ek(0 + firstequation_URB, 2 + firstequation_FRB + fDimension) += val_y;
        ek(2 + firstequation_URB, 0 + firstequation_FRB + fDimension) -= val_y;
        ek(1 + firstequation_URB, 0 + firstequation_FRB + fDimension) += val_z;
        ek(0 + firstequation_URB, 1 + firstequation_FRB + fDimension) -= val_z;
    }else{
        ek(0 + firstequation_URB, 0 + firstequation_FRB + fDimension) += val_y;
        ek(1 + firstequation_URB, 0 + firstequation_FRB + fDimension) -= val_x;
    }
    // Rotations x Rotations
    if(fDimension==3){
        ek(0 + firstequation_URB + fDimension, 0 + firstequation_FRB + fDimension) += val2_yz; 
        ek(1 + firstequation_URB + fDimension, 0 + firstequation_FRB + fDimension) -= val_xy; 
        ek(2 + firstequation_URB + fDimension, 0 + firstequation_FRB + fDimension) -= val_xz; 
        ek(0 + firstequation_URB + fDimension, 1 + firstequation_FRB + fDimension) -= val_xy; 
        ek(1 + firstequation_URB + fDimension, 1 + firstequation_FRB + fDimension) += val2_xz; 
        ek(2 + firstequation_URB + fDimension, 1 + firstequation_FRB + fDimension) -= val_yz; 
        ek(0 + firstequation_URB + fDimension, 2 + firstequation_FRB + fDimension) -= val_xz; 
        ek(1 + firstequation_URB + fDimension, 2 + firstequation_FRB + fDimension) -= val_yz; 
        ek(2 + firstequation_URB + fDimension, 2 + firstequation_FRB + fDimension) += val2_xy; 
    }else{
        ek(0 + firstequation_URB + fDimension, 0 + firstequation_FRB + fDimension) -= val2_xy;
        // std::cout<<"IntPoint: "<<delx<<"\n"; 
        // std::cout<<"val2_xy = "<<val2_xy<<"\n";
        // std::cout<<"weight = "<<weight<<"\n";
        // std::cout<<"ek(i,j) = "<<ek(0 + firstequation_URB + fDimension, 0 + firstequation_FRB + fDimension)<<std::endl;
    }
#ifdef PZDEBUG2
    ofstream out2("matrix2.txt");
    ek.Print("name",out2,EMathematicaInput);
    out2.flush();
    // if(!ek.IsSimetric()) DebugStop(); // TODO: delete this when done with debugging
#endif // PZDEBUG
}



void TPZMixedElasticityND::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
    if (datavec[0].fVecShapeIndex.size() == 0) {
        FillVecShapeIndex(datavec[0]);
    }
    int nspaces = datavec.size();
    switch(nspaces){
        case 3:
            Contribute_3spaces(datavec, weight, ek, ef);
            break;
        case 5:
        case 7:
            Contribute_5spaces(datavec, weight, ek, ef);
            break;
        default:
            DebugStop();
    }    
}


void TPZMixedElasticityND::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) {
    if (datavec[0].phi.Rows() != 0 && datavec[0].fShapeType != TPZMaterialData::EScalarShape) {
        DebugStop();
    }
    int ndisp = datavec[1].phi.Rows();

    const TPZVec<REAL> v_2 = bc.Val2();
    TPZFNMatrix<9, STATE> v_1 = bc.Val1();

    // Setting forcing function
    if (bc.HasForcingFunctionBC()) {
        TPZManVector<STATE, 3> res(fDimension);
        TPZFNMatrix<9, STATE> tens(fDimension, fDimension);
        bc.ForcingFunctionBC()(datavec[0].x, res, tens);
        v_2[0] = res[0];
        v_2[1] = res[1];
        if(fDimension == 3) v_2[2] = res[2];
    }

    // Setting the phis
    // E
    TPZFMatrix<REAL> &phiS = datavec[0].phi;

    int nshapeS;
    nshapeS = datavec[0].phi.Rows();


    //    TPZMaterialData::MShapeFunctionType shapetype = data.fShapeType;
    //    if(shapetype==data.EVecShape){
    //        ContributeVecShapeBC(data,weight,ek, ef,bc);
    //        return;
    //    }

    REAL R = datavec[0].x[0];
    if (R < 1.e-6) R = 1.e-6;
    
    int nstate = 2;
    if(fDimension == 3) nstate = 3;

    switch (bc.Type()) {
        case 0: // Dirichlet condition
        {
            for (int iq = 0; iq < nshapeS; iq++) {
                for (int idf = 0; idf < nstate; idf++) {
                    ef(nstate * iq + idf, 0) += v_2[idf] * phiS(iq, 0) * weight; // forced v2 displacement
                }
            }
        }
            break;

        case 1: // Neumann condition
        {
            for (int iq = 0; iq < nshapeS; iq++) {
                for (int jq = 0; jq < nshapeS; jq++) {
                    if (fAxisSymmetric) {
                        ek(2 * iq, 2 * jq) += TPZMaterial::fBigNumber * phiS(iq, 0) * phiS(jq, 0) * weight / R;
                        ek(2 * iq + 1, 2 * jq + 1) += TPZMaterial::fBigNumber * phiS(iq, 0) * phiS(jq, 0) * weight / R;
                    } else {
                        for(int idf = 0; idf < nstate; idf++)
                        {
                            ek(nstate * iq + idf, nstate * jq + idf) += TPZMaterial::fBigNumber * phiS(iq, 0) * phiS(jq, 0) * weight;
                        }
                    }
                }
                for(int idf = 0; idf < nstate; idf++)
                {
                    ef(nstate * iq + idf, 0) += TPZMaterial::fBigNumber * v_2[idf] * phiS(iq, 0) * weight; // normal stress in x direction
                }
            }
        }
            break;

        case 2: // Mixed condition
        {
            for (int iq = 0; iq < nshapeS; iq++) {
                for (int jq = 0; jq < nshapeS; jq++) {
                    if (fAxisSymmetric) {
                        ek(2 * iq, 2 * jq) += v_1(0, 0) * phiS(iq, 0) * phiS(jq, 0) * weight / R;
                        ek(2 * iq + 1, 2 * jq + 1) += v_1(1, 1) * phiS(iq, 0) * phiS(jq, 0) * weight / R;
                        ek(2 * iq + 1, 2 * jq) += v_1(1, 0) * phiS(iq, 0) * phiS(jq, 0) * weight / R;
                        ek(2 * iq, 2 * jq + 1) += v_1(0, 1) * phiS(iq, 0) * phiS(jq, 0) * weight / R;
                    } else {
                        for (int idf = 0; idf < nstate; idf++)
                        {
                            for (int jdf = 0; jdf < nstate; jdf++)
                            {
                                ek(nstate * iq + idf, nstate * jq + jdf) += v_1(idf, jdf) * phiS(iq, 0) * phiS(jq, 0) * weight;
                            }
                        }
                    }
                }
                for (int idf = 0; idf < nstate; idf++) {
                    ef(nstate * iq + idf, 0) += v_2[idf] * phiS(iq, 0) * weight; // normal stress
                }
            }
        }
            break;
        case 5:
            /// aplicando multiplicador de lagrange para impor deslocamento
            if (ndisp != 1) {
                DebugStop();
            }
            for (int idf = 0; idf < nstate; idf++)
            {
                if (v_1(idf, idf) > 1.e-1) {
                    ek(nstate * nshapeS + idf, nstate * nshapeS + idf) = 1.;
                }
            }
            break;
        default:
            DebugStop();
            // nulo introduzindo o BIGNUMBER pelos valores da condição
    } // 1 Val1 : a leitura é 00 01 10 11
    
}


/** Returns the variable index associated with the name. */
int TPZMixedElasticityND::VariableIndex(const std::string &name) const {
    if (!strcmp("displacement", name.c_str())) return 9;
    if (!strcmp("Displacement", name.c_str())) return 9; //function U ***
    if (!strcmp("DisplacementMem", name.c_str())) return 9;
    if (!strcmp("Pressure", name.c_str())) return 1;
    if (!strcmp("MaxStress", name.c_str())) return 2;
    if (!strcmp("PrincipalStress1", name.c_str())) return 3;
    if (!strcmp("PrincipalStress2", name.c_str())) return 4;
    if (!strcmp("SigmaX", name.c_str())) return 5;
    if (!strcmp("SigmaY", name.c_str())) return 6;
    if (!strcmp("SigmaZ", name.c_str())) return 12;
    if (!strcmp("TauXY", name.c_str())) return 8;
    if (!strcmp("TauXZ", name.c_str())) return 30;
    if (!strcmp("TauYZ", name.c_str())) return 31;
    if (!strcmp("Strain", name.c_str())) return 11; //Philippe
    if (!strcmp("SigmaZ", name.c_str())) return 12; //Philippe
    if (!strcmp("sig_x", name.c_str())) return 5;
    if (!strcmp("sig_y", name.c_str())) return 6;
    if (!strcmp("tau_xy", name.c_str())) return 8; //Cedric
    if (!strcmp("Displacement6", name.c_str())) return 7;
    if (!strcmp("Stress", name.c_str())) return 10; //function S ***
    if (!strcmp("Flux", name.c_str())) return 10;
    if (!strcmp("J2", name.c_str())) return 20;
    if (!strcmp("I1", name.c_str())) return 21;
    if (!strcmp("J2Stress", name.c_str())) return 20;
    if (!strcmp("I1Stress", name.c_str())) return 21;
    if (!strcmp("Alpha", name.c_str())) return 22;
    if (!strcmp("PlasticSqJ2", name.c_str())) return 22;
    if (!strcmp("PlasticSqJ2El", name.c_str())) return 22;
    if (!strcmp("YieldSurface", name.c_str())) return 27;
    if (!strcmp("NormalStress", name.c_str())) return 23;
    if (!strcmp("ShearStress", name.c_str())) return 24;
    if (!strcmp("NormalStrain", name.c_str())) return 25;
    if (!strcmp("ShearStrain", name.c_str())) return 26;
    if (!strcmp("Rotation", name.c_str())) return 27; //function P ***
    if(!strcmp("Young_Modulus",name.c_str()))        return 28;
    if(!strcmp("Poisson",name.c_str()))        return 29;
    if (!strcmp("ExactDisplacement", name.c_str())) return 33;
    if (!strcmp("ExactStress", name.c_str())) return 34;
    if (!strcmp("ElementSigmaError", name.c_str())) return 100;

    return TPZMaterial::VariableIndex(name);
}

/** Returns the number of variables associated with the variable indexed by var. */
int TPZMixedElasticityND::NSolutionVariables(int var) const {
    int nstate = fDimension;
    switch (var) {
        case 0:
            return nstate;
        case 1:
        case 2:
            return 1;
        case 3:
        case 4:
            return nstate;
        case 5:
        case 6:
        case 8:
            return 1;
        case 7:
            return 6;
        case 9:
            return 3;
        case 10: //Stress Tensor
            return nstate*nstate;
        case 11: //Strain Tensor
            return nstate*nstate;
            // SigZ
        case 12:
            return 1;
        case 20:
            return 1;
        case 21:
            return 1;
        case 22:
            return 1;
        case 23:
        case 24:
        case 25:
        case 26:
        case 27: // Rotation 
            return Dimension() == 3 ? 3 : 1;
        case 28:
        case 29:
        case 30:
        case 31:
            return 1;
        case 33:
        case 34:
            return nstate;
        case 100:
            return 1;
        default:
            return TPZMaterial::NSolutionVariables(var);
    }
}


/** @brief Returns the solution associated with the var index based on the finite element approximation */
void TPZMixedElasticityND::Solution(const TPZVec<TPZMaterialDataT<STATE>> &data, int var, TPZVec<STATE> &Solout) {
#ifdef PZDEBUG
    // if (data.size() != 3) {
    //     DebugStop();
    // }
#endif
    TPZManVector<REAL, 3> x = data[0].x;
    REAL R = x[0];
    if (R < 1.e-6) {
        R = 1.e-6;
    }
    TPZFNMatrix<9, STATE> sigma(3, 3, 0.), sigmah(3, 3, 0.), eps(3, 3, 0.);
    int dim = Dimension();
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < 3; j++) {
            if (fAxisSymmetric) {
                sigma(i, j) = data[0].sol[0][j + i * 3] / R;
            } else {
                sigma(i, j) = data[0].sol[0][j + i * 3];
            }
        }
    }
    TPZManVector<STATE, 3> disp(dim);
    if ( var == 9 && data[1].sol[0].size() != 0){
        for (int i = 0; i < dim; i++) {
            disp[i] = data[1].sol[0][i];
        }
    }
    
    TPZFNMatrix<9, STATE> antisym(dim, dim, 0.);
    antisym(0, 1) = data[2].sol[0][0];
    antisym(1, 0) = -antisym(0, 1);
    if(dim == 3)
    {
        antisym(0, 2) = data[2].sol[0][1];
        antisym(2, 0) = -antisym(0, 2);
        antisym(1, 2) = data[2].sol[0][2];
        antisym(2, 1) = -antisym(1, 2);
    }

    TElasticityAtPoint elast(fE_const,fnu_const);
    if(fElasticity)
    {
        //TPZManVector<REAL,3> result(2);
		TPZManVector<STATE, 3> result(dim);
        TPZFNMatrix<9,STATE> Dres(0,0);
        fElasticity(x, result, Dres);
        REAL E = result[0];
        REAL nu = result[1];
        TElasticityAtPoint modify(E,nu);
        elast = modify;
    }
    
    REAL mu = elast.fmu;
    REAL E = elast.fE;
    REAL nu = elast.fnu;
    
    if(var == 28)
    {
        Solout[0] = E;
        return ;
    }
    if(var == 29)
    {
        Solout[0] = elast.fnu;
        return;
    }
    

    REAL Pressure;

    //TPZManVector<REAL, 4> SIGMA(4, 0.), EPSZ(4, 0.);
	TPZManVector<STATE, 9> SIGMA(dim*dim, 0.), EPSZ(dim*dim, 0.);

    ToVoigt(sigma, SIGMA);

    ComputeDeformationVector(SIGMA, EPSZ, elast);

    FromVoigt(EPSZ, eps);

    if(dim ==2)
    {
        if (this->fPlaneStress == 1) {
            eps(2, 2) = -1. / E * nu * (sigma(0, 0) + sigma(1, 1));
        } else {
            sigma(2, 2) = nu * (sigma(0, 0) + sigma(1, 1));
            eps(2, 2) = 0;
        }
    }
    Pressure = -1/3. * (sigma(0, 0) + sigma(1, 1) + sigma(2,2));
    sigmah = sigma;
    sigmah(0, 0) += Pressure;
    sigmah(1, 1) += Pressure;
    sigmah(2, 2) += Pressure;
    // Displacement
    if (var == 9) {
        Solout[2] = 0.;
        for(int idf = 0; idf < dim; idf++)
        {
            Solout[idf] = disp[idf];
        }
        return;
    }
    // Pressure
    if (var == 1) {
        Solout[0] = Pressure;
        return;
    }
    // Sigmaz
    if (var == 12) {
        Solout[0] = sigma(2, 2);

        return;
    }
    // SigmaY
    if (var == 6) {
        Solout[0] = sigma(1, 1);
        return;
    }

    // SigmaX                
    if (var == 5) {
        Solout[0] = sigma(0, 0);
        return;
    }
    //  TauXY
    if (var == 8) {
        Solout[0] = 0.5 * (sigma(0, 1) + sigma(1, 0));
        return;
    }
    //  TauXZ
    if (var == 30) {
        Solout[0] = 0.5 * (sigma(0, 2) + sigma(2, 0));
        return;
    }
    //  TauYZ
    if (var == 31) {
        Solout[0] = 0.5 * (sigma(1, 2) + sigma(2, 1));
        return;
    }

    // Exact displacement                
    if (var == 33) {
        TPZVec<STATE> u_exact(fDimension,0.);
        TPZFMatrix<STATE> du_exact(fDimension,fDimension,0.);
        if (this->fExactSol) {
            this->fExactSol(data[0].x, u_exact, du_exact);
        }
        for (int idf = 0; idf < dim; idf++) {
            Solout[idf] = u_exact[idf];
        }
        return;
    }
    // Exact stress                
    if (var == 34) {
        
        TPZVec<STATE> u_exact(fDimension,0.);
        TPZFMatrix<STATE> du_exact(fDimension,fDimension,0.);
        if (this->fExactSol) {
            this->fExactSol(data[0].x, u_exact, du_exact);
        }
        // std::cout << "duexact = " << du_exact << std::endl;
        // For 2D only
        Solout[0] = du_exact(0,0);//Sigma x
        Solout[1] = du_exact(1,1);//Sigma y
        Solout[2] = (du_exact(1,0)+du_exact(0,1))/2;//Tau xy


        return;
    }




    //Strain
    if (var == 11) {
        Solout[Exx] = eps(0, 0);
        Solout[Eyx] = eps(1, 0);
        Solout[Exy] = eps(0, 1);
        Solout[Eyy] = eps(1, 1);
        if(dim == 3)
        {
            Solout[Exz] = eps(0, 2);
            Solout[Eyz] = eps(1, 2);
            Solout[Ezx] = eps(2, 0);
            Solout[Ezy] = eps(2, 1);
            Solout[Ezz] = eps(2, 2);
        }
        return;
    }

    //Stress
    if (var == 10) {
        Solout[Exx] = sigma(0, 0);
        Solout[Eyx] = sigma(1, 0);
        Solout[Exy] = sigma(0, 1);
        Solout[Eyy] = sigma(1, 1);
        if(dim == 3)
        {
            Solout[Exz] = sigma(0, 2);
            Solout[Eyz] = sigma(1, 2);
            Solout[Ezx] = sigma(2, 0);
            Solout[Ezy] = sigma(2, 1);
            Solout[Ezz] = sigma(2, 2);
        }
        return;
    }

    //I1
    if (var == 21) {
        Solout[0] = sigma(0, 0) + sigma(1, 1) + sigma(2, 2);
        return;
    }
    //J2
    if (var == 20) {
        Solout[0] = sigmah(1, 1) * sigmah(0, 0) + sigmah(1, 1) * sigma(2, 2)
         + sigmah(0,0) * sigmah(2,2) - sigmah(1, 0) * sigmah(1, 0) - sigmah(0, 1) * sigmah(0, 1);
        return;
    }
    //NormalStrain?

    //PrincipalStrain1
    if (var == 3) {
        Solout[0] = Pressure + sqrt(0.25 * (sigma(0, 0) - sigma(1, 1))*(sigma(0, 0) - sigma(1, 1)) + sigma(1, 2));
        return;
    }
    //Rotation
    if (var == 27) {
        Solout[0] = antisym(0, 1);
        if(dim == 3)
        {
            Solout[1] = antisym(0, 2);
            Solout[2] = antisym(1, 2);
        }
        return;
    }
}


////////////////////////////////////////////////////////////////////

STATE TPZMixedElasticityND::Inner(TPZFMatrix<STATE> &S, TPZFMatrix<STATE> &T) {
    //inner product of two tensors
#ifdef DEBUG
    if (S.Rows() != S.Cols() || T.Cols() != T.Rows() || S.Rows() != T.Rows()) {
        DebugStop();
    }
#endif

    STATE Val = 0;
    int sr = S.Rows();
    for (int i = 0; i < sr; i++) {
        for (int j = 0; j < sr; j++) {
            Val += S(i, j) * T(i, j);
        }
    }
    return Val;
}

////////////////////////////////////////////////////////////////////

template <typename TVar>
TVar TPZMixedElasticityND::InnerVec(const TPZVec<TVar> &S, const TPZVec<TVar> &T) {
    //inner product of two vectors
#ifdef DEBUG
    if (S.size() != T.size()) {
        DebugStop();
    }
#endif
    TVar Val = 0;
    for (int i = 0; i < S.size(); i++) {
        Val += S[i] * T[i];
    }
    return Val;
}

////////////////////////////////////////////////////////////////////
STATE TPZMixedElasticityND::Tr(TPZFMatrix<REAL> &GradU) {
#ifdef DEBUG
    if (GradU.Rows() != GradU.Cols()) {
        DebugStop();
    }
#endif

    STATE Val = 0.;

    int grr = GradU.Rows();
    for (int i = 0; i < grr; i++) {
        Val += GradU(i, i);
    }

    return Val;
}


/// transform a H1 data structure to a vector data structure

void TPZMixedElasticityND::FillVecShapeIndex(TPZMaterialData &data) {
    data.fDeformedDirections.Resize(fDimension, fDimension);
    data.fDeformedDirections.Identity();
    data.fVecShapeIndex.Resize(fDimension * data.phi.Rows());
    for (int d = 0; d < fDimension; d++) {
        for (int i = 0; i < data.phi.Rows(); i++) {
            data.fVecShapeIndex[i * fDimension + d].first = d;
            data.fVecShapeIndex[i * fDimension + d].second = i;
        }
    }
}


int TPZMixedElasticityND::NEvalErrors() const {
    return 7;
}

void TPZMixedElasticityND::Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors) {
    //values[0] = 0.;
    //TPZManVector<REAL, 4> SigmaV(4, 0.), sigma_exactV(4, 0.), eps_exactV(4, 0.), EPSZV(4, 0.);
    TPZVec<STATE> u_exact(fDimension,0.);
    TPZFMatrix<STATE> du_exact(fDimension,fDimension,0.);
    TPZManVector<STATE> divsigma(fDimension,0.);
    if (this->fExactSol) {
        this->fExactSol(data[0].x, u_exact, du_exact);
    }
    if (this->fForcingFunction) {
        this->fForcingFunction(data[0].x, divsigma);
    }

    int nstate = fDimension;
    int matdim = nstate*nstate;
	TPZManVector<STATE, 9> SigmaV(matdim, 0.), sigma_exactV(matdim, 0.), eps_exactV(matdim, 0.), EPSZV(matdim, 0.);
	TPZFNMatrix<9, STATE> sigma(nstate, nstate, 0.), eps(nstate, nstate, 0.), grad(nstate, nstate, 0.);
    TPZFNMatrix<9, STATE> eps_exact(nstate, nstate, 0.);
    TPZManVector<REAL, 3> x = data[0].x;
    int dim = Dimension();
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            sigma(i, j) = data[0].sol[0][j + i * 3];
            if (fAxisSymmetric) {
                sigma(i, j) /= x[0];
            }
        }
    }
    ToVoigt(sigma, SigmaV);
    TPZManVector<STATE, 3> divSigma(dim,0.);
    for (int i = 0; i < dim; i++) {
        divSigma[i] = data[0].divsol[0][i];
    }

    TPZManVector<STATE, 3> disp(dim);
    for (int i = 0; i < dim; i++) {
        disp[i] = data[1].sol[0][i];
    }

#ifdef LOG4CXX
    if (logdata->isDebugEnabled())
 {
        std::stringstream sout;
        sout << "DISP*************************** = " << disp << std::endl;
        sigma.Print("sigma************************************ = ", sout, EMathematicaInput);
        LOGPZ_DEBUG(logdata, sout.str())
    }
#endif


    eps_exact(0, 0) = du_exact(0, 0);
    eps_exact(1, 0) = eps_exact(0, 1) = 0.5 * (du_exact(0, 1) + du_exact(1, 0));
    eps_exact(1, 1) = du_exact(1, 1);
    if(dim == 3)
    {
        eps_exact(2, 2) = du_exact(2, 2);
        eps_exact(2, 0) = eps_exact(0, 2) = 0.5 * (du_exact(0, 2) + du_exact(2, 0));
        eps_exact(1, 2) = eps_exact(2, 1) = 0.5 * (du_exact(2, 1) + du_exact(1, 2));
    }
    ToVoigt(eps_exact, eps_exactV);
    
    TElasticityAtPoint elast(fE_const,fnu_const);
    if(fElasticity)
    {
        //TPZManVector<REAL,3> result(2);
		TPZManVector<STATE, 3> result(2);
        TPZFNMatrix<4,STATE> Dres(0,0);
        fElasticity(x, result, Dres);
        REAL E = result[0];
        REAL nu = result[1];
        TElasticityAtPoint modify(E,nu);
        elast = modify;
    }
    
    ComputeStressVector(eps_exactV, sigma_exactV, elast);

    TPZManVector<STATE,3> rotation(1,0.), rotationExact(1,0.);
    if(dim == 3) {
        rotation.Resize(fDimension,0.);
        rotationExact.Resize(fDimension,0.);
    }
    rotation[0] = data[2].sol[0][0];
    rotationExact[0] = 0.5*(du_exact(1, 0)-du_exact(0, 1));
    if(dim == 3)
    {
        rotation[1] = data[2].sol[0][1];
        rotationExact[1] = 0.5*(du_exact(2, 0)-du_exact(0, 2));
        rotation[2] = data[2].sol[0][2];
        rotationExact[2] = 0.5*(du_exact(2, 1)-du_exact(1, 2));

    }
    // std::cout << "Rotation = "<< rotation << std::endl;
    // std::cout << "rotationExact = "<< rotationExact << std::endl << std::endl;

#ifdef PZDEBUG
    {
        TPZManVector<STATE, 9> eps_again(matdim);
        ComputeDeformationVector(sigma_exactV, eps_again, elast);
        for (int i = 0; i < matdim; i++) {
            if (abs(eps_again[i] - eps_exactV[i]) > 1.e-8) {
                DebugStop();
            }
        }
    }
#endif
    
    //L2_Error: for the displacement (u)
    errors[3] = 0.;
    for (int idf = 0; idf < dim; idf++) {
        errors[3] += (disp[idf] - u_exact[idf])*(disp[idf] - u_exact[idf]);
    }
    //Energe norm
    //TPZManVector<REAL,4> SIGMA(4,0.) , EPSZ(4,0.);

    
    ToVoigt(sigma, SigmaV);

    ComputeDeformationVector(SigmaV, EPSZV, elast);

    errors[0] = 0.;
    errors[1] = 0.;
    for (int i = 0; i < matdim; i++) {
        //L2_Error: for the stress tensor (sigma)
        errors[0] += (SigmaV[i] - sigma_exactV[i])*(SigmaV[i] - sigma_exactV[i]);
        
        //Energy_Error: for the stress tensor (sigma)
        errors[1] += (SigmaV[i] - sigma_exactV[i])*(EPSZV[i] - eps_exactV[i]);
    }
    if (errors[1] < 0.) {
        std::cout << "I should stop \n";
    }

    TPZManVector<STATE,3> divSigmaExact(fDimension,0.);
    if(HasForcingFunction())
    {
        fForcingFunction(data[0].x, divSigmaExact);
        for(int i=0; i<fDimension; i++) divSigmaExact[i] *= -1.;
    }
    errors[2] = 0.;
    errors[4] = 0.;
    for (int idf = 0; idf < fDimension; idf++)
    {
        //L2_Error: divergent of the stress tensor (div(sigma))
        
        errors[2] += (divSigma[idf]-divSigmaExact[idf])*(divSigma[idf]-divSigmaExact[idf]);
        
    }
    int nrot = 1;
    if(fDimension == 3) nrot = 3;
    for(int ir = 0; ir<nrot; ir++)
    {
        //L2_Error: for the rotation (q)
        errors[4] += (rotation[ir]-rotationExact[ir])*(rotation[ir]-rotationExact[ir]);
    }
    //L2_Error: for the symmetry measure (asym)
    errors[5] = pow(SigmaV[Exy]-SigmaV[Eyx],2);
    if(fDimension == 3)
    {
        errors[5] += pow(SigmaV[Exz]-SigmaV[Ezx],2);
        errors[5] += pow(SigmaV[Ezy]-SigmaV[Eyz],2);

    }
    
    //Energy_Norm: for the exact displacement solution
    errors[6] = 0.;
    for(int i=0; i<matdim; i++)
    {
        errors[6] += eps_exactV[i]*sigma_exactV[i];
    }
    //    std::cout << "x " << data[0].x << std::endl;
    //    std::cout << "disp " << u_exact << std::endl;
    //    du_exact.Print("du ",std::cout);
    //    std::cout << "sigma_exact " << sigma_exactV << std::endl;
    //    std::cout << errors << std::endl;
    //or we can compute as this methods
    //TPZFMatrix<STATE> MatrixElast(4,4,0.);
    //ElasticityModulusTensor(MatrixElast);
    //errors[1] = 0;
    //for(int i==0; i<4； i++)
    //{
    //for(int j=0; j<4; j++)
    //{
    //errors[1]+=SIGMA[j]*MatrixElast(j,i);
    //}
    //}

    // SemiH1 norm
    //TPZFNMatrix<4,REAL> antisym(2,2,0.);
    //antisym(0,1) = data[2].sol[0][0];
    //antisym(1,0) = -antisym(0,1);
    //TPZManVector<REAL,4> P(4,0.),DU(4,0.);
    //ToVoigt(antisym,P);
    //ToVoigt(du_exact,DU);
    //errors[2]=0;
    //for(int i=0;i<4;i++)
    //{
    //   errors[2]+=(P[i]+EPSZ[i]-DU[i])*(P[i]+EPSZ[i]-DU[i]);
    //}
    //values[0] = calculo do erro estimado em norma Energia
    //values[0] = fE*(sigx*sigx + sigy*sigy + 2*fnu*sigx*sigy)/(1-fnu*fnu);
    //values[0] = (values[0] + .5*fE*sigxy*sigxy/(1+fnu));

    //values[1] : erro em norma L2 em tens�s
    //values[1] = sigx*sigx + sigy*sigy + sigxy*sigxy;

    //values[1] : erro em norma L2 em deslocamentos
    //values[1] = pow((REAL)fabs(u[0] - u_exact[0]),(REAL)2.0)+pow((REAL)fabs(u[1] - u_exact[1]),(REAL)2.0);

    //values[2] : erro estimado na norma H1
    //REAL SemiH1 =0.;
    //for(int i = 0; i < 2; i++) for(int j = 0; j < 2; j++) SemiH1 += (du(i,j) - du_exact(i,j)) * (du(i,j) - du_exact(i,j));
    //values[2] = values[1] + SemiH1;
}


TPZMixedElasticityND::TPZMixedElasticityND(const TPZMixedElasticityND &copy) :
TBase(copy),
fE_const(copy.fE_const),
fnu_const(copy.fnu_const),
flambda_const(copy.flambda_const),
fmu_const(copy.fmu_const),
fForce(copy.fForce),
fPlaneStress(copy.fPlaneStress),
fDimension(copy.fDimension),
fMatrixA(copy.fMatrixA),
fElasticity(copy.fElasticity),
fAxisSymmetric(copy.fAxisSymmetric){
}

int TPZMixedElasticityND::ClassId() const {
    return Hash("TPZMixedElasticityND") ^ TBase::ClassId() << 1;
}

void TPZMixedElasticityND::Read(TPZStream &buf, void *context) {
    TPZMaterial::Read(buf, context);
    buf.Read(&fE_const, 1);
    buf.Read(&fnu_const, 1);
    fForce.Resize(2, 0.);
    buf.Read(&fForce[0], 2);
    buf.Read(&fPlaneStress, 1);
}

void TPZMixedElasticityND::Write(TPZStream &buf, int withclassid) const {
    TPZMaterial::Write(buf, withclassid);
    buf.Write(&fE_const, 1);
    buf.Write(&fnu_const, 1);

    buf.Write(&fForce[0], 2);
    buf.Write(&fPlaneStress, 1);
}

