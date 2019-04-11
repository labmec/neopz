#include "TPZSolveMatrix.h"
#include "TPZTensor.h"
#include "pzmatrix.h"
#include <stdlib.h>
#include "TPZTensor.h"
#include "TPZVTKGeoMesh.h"
#include "pzintel.h"

#ifdef USING_MKL
#include <mkl.h>
#include <algorithm>
#endif

//Spectral decomposition
void TPZSolveMatrix::Multiplicity1(double *sigma, double eigenvalue, double *eigenvector) {
    TPZVec<REAL> det(3);
    det[0] = (sigma[0] - eigenvalue)*(sigma[1] - eigenvalue) - sigma[3]*sigma[3];
    det[1] = (sigma[0] - eigenvalue)*(sigma[2] - eigenvalue);
    det[2] = (sigma[1] - eigenvalue)*(sigma[2] - eigenvalue);

    REAL maxdet = fabs(det[0]);
    for (int i = 1; i < 3; i++) {
        if (fabs(det[i]) > fabs(maxdet)) {
            maxdet = fabs(det[i]);
        }
    }
    TPZVec<REAL> v(3);
    if (maxdet == fabs(det[0])) {
        v[0] = 0;
        v[1] = 0;
        v[2] = 1;

    }
    else if (maxdet == fabs(det[1])) {
        v[0] = 1/det[1]*(-(sigma[2] - eigenvalue)*sigma[3]);
        v[1] = 1;
        v[2] = 0;

    }
    else {
        v[0] = 1;
        v[1] = 1/det[2]*(-(sigma[2] - eigenvalue)*sigma[3]);
        v[2] = 0;
    }
    REAL norm = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    eigenvector[0] = v[0]/norm;
    eigenvector[1] = v[1]/norm;
    eigenvector[2] = v[2]/norm;
}

void TPZSolveMatrix::Multiplicity2(double *sigma, double eigenvalue, double *eigenvector1, double *eigenvector2) {
    TPZVec<REAL> x(3);
    x[0] = sigma[0] - eigenvalue;
    x[1] = sigma[1] - eigenvalue;
    x[2] = sigma[2] - eigenvalue;

    REAL maxx = fabs(x[0]);
    for (int i = 1; i < 3; i++) {
        if (fabs(x[i]) > fabs(maxx)) {
            maxx = fabs(x[i]);
        }
    }

    TPZVec<REAL> v1(3);
    TPZVec<REAL> v2(3);

    if (maxx == fabs(x[0])) {
        v1[0] = -sigma[3]/x[0];
        v1[1] = 1;
        v1[2] = 0;

        v2[0] = 0;
        v2[1] = 0;
        v2[2] = 1;

    }
    else if (maxx == fabs(x[1])) {
        v1[0] = 1;
        v1[1] = -sigma[3]/x[1];
        v1[2] = 0;

        v2[0] = 0;
        v2[1] = 0;
        v2[2] = 1;

    }
    else {
        v1[0] = 1;
        v1[1] = 0;
        v1[2] = 0;

        v2[0] = 0;
        v2[1] = 1;
        v2[2] = 0;

    }
    REAL norm1 = sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);
    REAL norm2 = sqrt(v2[0]*v2[0] + v2[1]*v1[1] + v2[2]*v2[2]);

    eigenvector1[0] = v1[0]/norm1;
    eigenvector1[1] = v1[1]/norm1;
    eigenvector1[2] = v1[2]/norm1;

    eigenvector2[0] = v2[0]/norm2;
    eigenvector2[1] = v2[1]/norm2;
    eigenvector2[2] = v2[2]/norm2;
}

void TPZSolveMatrix::Eigenvectors(double *sigma, double *eigenvalues, double *eigenvectors, double &maxel) {
    sigma[0]*=maxel;
    sigma[1]*=maxel;
    sigma[2]*=maxel;
    sigma[3]*=maxel;

    if ((eigenvalues[0] == eigenvalues[1]) && (eigenvalues[1] == eigenvalues[2])) {
        eigenvectors[0] = 1.;
        eigenvectors[1] = 0.;
        eigenvectors[2] = 0.;

        eigenvectors[3] = 0.;
        eigenvectors[4] = 1.;
        eigenvectors[5] = 0.;

        eigenvectors[6] = 0.;
        eigenvectors[7] = 0.;
        eigenvectors[8] = 1.;
    }
    else {
        if (eigenvalues[0] != eigenvalues[1] && eigenvalues[0] != eigenvalues[2]) {
            Multiplicity1(sigma, eigenvalues[0], &eigenvectors[0]);
        } else if (eigenvalues[0] == eigenvalues[1]) {
            Multiplicity2(sigma, eigenvalues[0], &eigenvectors[0], &eigenvectors[3]);
        } else if (eigenvalues[0] == eigenvalues[2]) {
            Multiplicity2(sigma, eigenvalues[0], &eigenvectors[0], &eigenvectors[6]);
        }
        if (eigenvalues[1] != eigenvalues[0] && eigenvalues[1] != eigenvalues[2]) {
            Multiplicity1(sigma, eigenvalues[1], &eigenvectors[3]);
        } else if (eigenvalues[1] == eigenvalues[2]) {
            Multiplicity2(sigma, eigenvalues[1], &eigenvectors[3], &eigenvectors[6]);
        }
        if (eigenvalues[2] != eigenvalues[0] && eigenvalues[2] != eigenvalues[1]) {
            Multiplicity1(sigma, eigenvalues[2], &eigenvectors[6]);
        }
    }
}

void TPZSolveMatrix::Normalize(double *sigma, double &maxel) {
    maxel = sigma[0];
    for (int i = 1; i < 4; i++) {
        if (fabs(sigma[i]) > fabs(maxel)) {
            maxel = sigma[i];
        }
    }
    for (int i = 0; i < 4; i++) {
        sigma[i] /= maxel;
    }
}

void TPZSolveMatrix::Interval(double *sigma, double *interval) {
    TPZVec<REAL> lower_vec(3);
    TPZVec<REAL> upper_vec(3);

    //row 1 |sigma_xx sigma_xy 0|
    lower_vec[0] = sigma[0] - fabs(sigma[3]);
    upper_vec[0] = sigma[0] + fabs(sigma[3]);

    //row 2 |sigma_xy sigma_yy 0|
    lower_vec[1] = sigma[1] - fabs(sigma[3]);
    upper_vec[1] = sigma[1] + fabs(sigma[3]);

    //row 3 |0 0 sigma_zz|
    lower_vec[2] = sigma[2];
    upper_vec[2] = sigma[2];

    interval[0] = upper_vec[0];
    interval[1] = lower_vec[0];

    for (int i = 1; i < 3; i++) {
        if (upper_vec[i] > interval[0]) { //upper interval
            interval[0] = upper_vec[i];
        }

        if (lower_vec[i] < interval[1]) { //lower interval
            interval[1] = lower_vec[i];
        }
    }
}

void TPZSolveMatrix::NewtonIterations(double *interval, double *sigma, double *eigenvalues, double &maxel) {
    int numiterations = 20;
    REAL tol = 10e-12;

    REAL res, f, df, x;
    int it;

    for (int i = 0; i < 2; i++) {
        x = interval[i];
        it = 0;

        f = sigma[0] * sigma[1] - x * (sigma[0] + sigma[1]) + x * x - sigma[3] * sigma[3];
        res = abs(f);

        while (it < numiterations && res > tol) {
            df = -sigma[0] - sigma[1] + 2 * x;

            x -= f / df;
            f = sigma[0] * sigma[1] - x * (sigma[0] + sigma[1]) + x * x - sigma[3] * sigma[3];
            res = abs(f);
            it++;
        }
        eigenvalues[i] = x;

    }
    eigenvalues[2] = sigma[0] + sigma[1] + sigma[2] - eigenvalues[0] - eigenvalues[1];

    eigenvalues[0] *= maxel;
    eigenvalues[1] *= maxel;
    eigenvalues[2] *= maxel;

//    std::sort(eigenvalues, eigenvalues+3, [](int i, int j) { return i > j; }); //store eigenvalues in descending order (absolute value)
    std::sort(eigenvalues, eigenvalues+3, greater<REAL>()); //store eigenvalues in descending order (absolute value)

}

//Project Sigma
bool TPZSolveMatrix::PhiPlane(double *eigenvalues, double *sigma_projected) {
    REAL mc_phi = fMaterialData.FrictionAngle();
    REAL mc_cohesion = fMaterialData.Cohesion();

    const REAL sinphi = sin(mc_phi);
    const REAL cosphi = cos(mc_phi);

    REAL phi = eigenvalues[0] - eigenvalues[2] + (eigenvalues[0] + eigenvalues[2]) * sinphi - 2. * mc_cohesion *cosphi;

    sigma_projected[0] = eigenvalues[0];
    sigma_projected[1] = eigenvalues[1];
    sigma_projected[2] = eigenvalues[2];

    bool check_validity = (fabs(phi) < 1.e-12) || (phi < 0.0);
    return check_validity;
}

bool TPZSolveMatrix::ReturnMappingMainPlane(double *eigenvalues, double *sigma_projected, double &m_hardening) {
    REAL mc_phi = fMaterialData.FrictionAngle();
    REAL mc_psi = mc_phi;
    REAL mc_cohesion = fMaterialData.Cohesion();
    REAL G = fMaterialData.ElasticResponse().G();
    REAL K = fMaterialData.ElasticResponse().K();

    const REAL sinphi = sin(mc_phi);
    const REAL sinpsi = sin(mc_psi);
    const REAL cosphi = cos(mc_phi);
    const REAL sinphi2 = sinphi*sinphi;
    const REAL cosphi2 = 1. - sinphi2;
    const REAL constA = 4. * G *(1. + sinphi * sinpsi / 3.) + 4. * K * sinphi*sinpsi;

    REAL phi = eigenvalues[0] - eigenvalues[2]+(eigenvalues[0] + eigenvalues[2]) * sinphi - 2. * mc_cohesion*cosphi;

    REAL gamma = 0;
    int n_iterations = 30;
    for (int i = 0; i < n_iterations; i++) {
        double jac = -constA - 4. * cosphi2 * 0; // H=0
        double delta_gamma = - phi / jac;
        gamma += delta_gamma;
        phi = eigenvalues[0] - eigenvalues[2]+(eigenvalues[0] + eigenvalues[2]) * sinphi - 2. * mc_cohesion * cosphi - constA * gamma;
        if (fabs(phi) < 1.e-12) {
            break;
        }
    }

    sigma_projected[0] -= (2. * G *(1 + sinpsi / 3.) + 2. * K * sinpsi) * gamma;
    sigma_projected[1] += (4. * G / 3. - K * 2.) * sinpsi * gamma;
    sigma_projected[2] += (2. * G * (1 - sinpsi / 3.) - 2. * K * sinpsi) * gamma;

    m_hardening += gamma * 2. * cosphi;

    bool check_validity = (eigenvalues[0] > eigenvalues[1] || fabs(eigenvalues[0]-eigenvalues[1]) < 1.e-12) && (eigenvalues[1] > eigenvalues[2] || fabs(eigenvalues[1]-eigenvalues[2]) < 1.e-12);
    return check_validity;
}

bool TPZSolveMatrix::ReturnMappingRightEdge(double *eigenvalues, double *sigma_projected, double &m_hardening) {
    REAL mc_phi = fMaterialData.FrictionAngle();
    REAL mc_psi = mc_phi;
    REAL mc_cohesion = fMaterialData.Cohesion();
    REAL G = fMaterialData.ElasticResponse().G();
    REAL K = fMaterialData.ElasticResponse().K();

    const REAL sinphi = sin(mc_phi);
    const REAL sinpsi = sin(mc_psi);
    const REAL cosphi = cos(mc_phi);

    TPZVec<REAL> gamma(2, 0.), phi(2, 0.), sigma_bar(2, 0.), ab(2, 0.);

    TPZVec<TPZVec<REAL>> jac(2), jac_inv(2);
    for (int i = 0; i < 2; i++) {
        jac[i].Resize(2, 0.);
        jac_inv[i].Resize(2, 0.);
    }

    sigma_bar[0] = eigenvalues[0] - eigenvalues[2]+(eigenvalues[0] + eigenvalues[2]) * sinphi;
    sigma_bar[1] = eigenvalues[0] - eigenvalues[1] + (eigenvalues[0] + eigenvalues[1]) * sinphi;

    phi[0] = sigma_bar[0] - 2. * cosphi * mc_cohesion;
    phi[1] = sigma_bar[1] - 2. * cosphi * mc_cohesion;

    ab[0] = 4. * G * (1 + sinphi * sinpsi / 3.) + 4. * K * sinphi * sinpsi;
    ab[1] = 2. * G * (1. + sinphi + sinpsi - sinphi * sinpsi / 3.) + 4. * K * sinphi * sinpsi;

    int n_iterations = 30;
    for (int i = 0; i < n_iterations; i++) {

        jac[0][0] = -ab[0];
        jac[1][0] = -ab[1];
        jac[0][1] = -ab[1];
        jac[1][1] = -ab[0];

        double det_jac = jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0];

        jac_inv[0][0] = jac[1][1] / det_jac;
        jac_inv[1][0] = -jac[1][0] / det_jac;
        jac_inv[0][1] = -jac[0][1] / det_jac;
        jac_inv[1][1] = jac[0][0] / det_jac;

        gamma[0] -= (jac_inv[0][0] * phi[0] + jac_inv[0][1] * phi[1]);
        gamma[1] -= (jac_inv[1][0] * phi[0] + jac_inv[1][1] * phi[1]);

        phi[0] = sigma_bar[0] - ab[0] * gamma[0] - ab[1] * gamma[1] - 2. * cosphi * mc_cohesion;
        phi[1] = sigma_bar[1] - ab[1] * gamma[0] - ab[0] * gamma[1] - 2. * cosphi * mc_cohesion;

        double res = (fabs(phi[0]) + fabs(phi[1]));

        if (fabs(res) < 1.e-12) {
            break;
        }
    }

    sigma_projected[0] -= (2. * G * (1 + sinpsi / 3.) + 2. * K * sinpsi) * (gamma[0] + gamma[1]);
    sigma_projected[1] += ((4. * G / 3. - K * 2.) * sinpsi) * gamma[0] + (2. * G * (1. - sinpsi / 3.) - 2. * K * sinpsi) * gamma[1];
    sigma_projected[2] += (2. * G * (1 - sinpsi / 3.) - 2. * K * sinpsi) * gamma[0] + ((4. * G / 3. - 2. * K) * sinpsi) * gamma[1];

    m_hardening += (gamma[0] + gamma[1]) * 2. * cosphi;

    bool check_validity = (eigenvalues[0] > eigenvalues[1] || fabs(eigenvalues[0]-eigenvalues[1]) < 1.e-12) && (eigenvalues[1] > eigenvalues[2] || fabs(eigenvalues[1]-eigenvalues[2]) < 1.e-12);
    return check_validity;
}

bool TPZSolveMatrix::ReturnMappingLeftEdge(double *eigenvalues, double *sigma_projected, double &m_hardening) {
    REAL mc_phi = fMaterialData.FrictionAngle();
    REAL mc_psi = mc_phi;
    REAL mc_cohesion = fMaterialData.Cohesion();
    REAL G = fMaterialData.ElasticResponse().G();
    REAL K = fMaterialData.ElasticResponse().K();

    const REAL sinphi = sin(mc_phi);
    const REAL sinpsi = sin(mc_psi);
    const REAL cosphi = cos(mc_phi);
    const REAL sinphi2 = sinphi*sinphi;
    const REAL cosphi2 = 1. - sinphi2;

    TPZVec<REAL> gamma(2, 0.), phi(2, 0.), sigma_bar(2, 0.), ab(2, 0.);

    TPZVec<TPZVec<REAL>> jac(2), jac_inv(2);
    for (int i = 0; i < 2; i++) {
        jac[i].Resize(2, 0.);
        jac_inv[i].Resize(2, 0.);
    }

    sigma_bar[0] = eigenvalues[0] - eigenvalues[2]+(eigenvalues[0] + eigenvalues[2]) * sinphi;
    sigma_bar[1] = eigenvalues[1] - eigenvalues[2]+(eigenvalues[1] + eigenvalues[2]) * sinphi;

    ab[0] = 4. * G * (1 + sinphi * sinpsi / 3.) + 4. * K * sinphi * sinpsi;
    ab[1] = 2. * G * (1. - sinphi - sinpsi - sinphi * sinpsi / 3.) + 4. * K * sinphi * sinpsi;

    phi[0] = sigma_bar[0] - 2. * cosphi * mc_cohesion;
    phi[1] = sigma_bar[1] - 2. * cosphi * mc_cohesion;

    int n_iterations = 30;
    for (int i = 0; i < n_iterations; i++) {

        jac[0][0] = -ab[0] - 4. * cosphi2 * 0;
        jac[1][0] = -ab[1] - 4. * cosphi2 * 0;
        jac[0][1] = -ab[1] - 4. * cosphi2 * 0;
        jac[1][1] = -ab[0] - 4. * cosphi2 * 0;

        REAL det_jac = jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0];

        jac_inv[0][0] = jac[1][1] / det_jac;
        jac_inv[1][0] = -jac[1][0] / det_jac;
        jac_inv[0][1] = -jac[0][1] / det_jac;
        jac_inv[1][1] = jac[0][0] / det_jac;

        gamma[0] -= (jac_inv[0][0] * phi[0] + jac_inv[0][1] * phi[1]);
        gamma[1] -= (jac_inv[1][0] * phi[0] + jac_inv[1][1] * phi[1]);

        phi[0] = sigma_bar[0] - ab[0] * gamma[0] - ab[1] * gamma[1] - 2. * cosphi * mc_cohesion;
        phi[1] = sigma_bar[1] - ab[1] * gamma[0] - ab[0] * gamma[1] - 2. * cosphi * mc_cohesion;

        REAL res = (fabs(phi[0]) + fabs(phi[1]));

        if (fabs(res) < 1.e-12) {
            break;
        }
    }

    sigma_projected[0] += -(2. * G * (1 + sinpsi / 3.) + 2. * K * sinpsi) * gamma[0] + ((4. * G / 3. - 2. * K) * sinpsi) * gamma[1];
    sigma_projected[1] += ((4. * G / 3. - K * 2.) * sinpsi) * gamma[0] - (2. * G * (1. + sinpsi / 3.) + 2. * K * sinpsi) * gamma[1];
    sigma_projected[2] += (2. * G * (1 - sinpsi / 3.) - 2. * K * sinpsi) * (gamma[0] + gamma[1]);

    m_hardening += (gamma[0] + gamma[1]) * 2. * cosphi;

    bool check_validity = (eigenvalues[0] > eigenvalues[1] || fabs(eigenvalues[0]-eigenvalues[1]) < 1.e-12) && (eigenvalues[1] > eigenvalues[2] || fabs(eigenvalues[1]-eigenvalues[2]) < 1.e-12);
    return check_validity;
}

void TPZSolveMatrix::ReturnMappingApex(double *eigenvalues, double *sigma_projected, double &m_hardening) {
    REAL mc_phi = fMaterialData.FrictionAngle();
    REAL mc_psi = mc_phi;
    REAL mc_cohesion = fMaterialData.Cohesion();
    REAL K = fMaterialData.ElasticResponse().K();

    const REAL cotphi = 1. / tan(mc_phi);

    REAL ptrnp1 = 0.;
    for (int i = 0; i < 3; i++) {
        ptrnp1 += eigenvalues[i];
    }
    ptrnp1 /= 3.;

    REAL DEpsPV = 0.;
    REAL alpha = cos(mc_phi) / sin(mc_psi);
    REAL res = mc_cohesion * cotphi - ptrnp1;
    REAL pnp1;

    int n_iterations = 30;
    for (int i = 0; i < n_iterations; i++) {
        const REAL jac = K; //H = 0
        DEpsPV -= res / jac;

        pnp1 = ptrnp1 - K * DEpsPV;
        res = mc_cohesion * cotphi - pnp1;

        if (fabs(res) < 1.e-12) {
            break;
        }
    }

    m_hardening += alpha * DEpsPV;
    for (int i = 0; i < 3; i++) {
        sigma_projected[i] = pnp1;
    }
}

//Gather solution
void TPZSolveMatrix::GatherSolution(TPZFMatrix<REAL> &global_solution, TPZFMatrix<REAL> &gather_solution) {
    int64_t n_globalsol = fCmesh->Dimension()*fNphis;

    gather_solution.Resize(n_globalsol,1);
    gather_solution.Zero();

    cblas_dgthr(n_globalsol, global_solution, &gather_solution(0,0), &fIndexes[0]);

}

//Strain
void TPZSolveMatrix::DeltaStrain(TPZFMatrix<REAL> &expandsolution, TPZFMatrix<REAL> &delta_strain) {
    int64_t nelem = fRowSizes.size();
    int64_t n_globalsol = fCmesh->Dimension()*fNphis;

    delta_strain.Resize(2*n_globalsol,1);
    delta_strain.Zero();

    for (int64_t iel = 0; iel < nelem; iel++) {
        for (int i = 0; i < fRowSizes[iel]; i++) {
                for (int k = 0; k < fColSizes[iel]; k++) {
                    delta_strain(i + fRowFirstIndex[iel], 0) += fStorage[k * fRowSizes[iel] + i + fMatrixPosition[iel]] * expandsolution(k + fColFirstIndex[iel],0);
                    delta_strain(i + fRowFirstIndex[iel] + n_globalsol, 0) += fStorage[k * fRowSizes[iel] + i + fMatrixPosition[iel]] * expandsolution(k + fColFirstIndex[iel] + n_globalsol/2,0);
                }
        }
    }
}

void TPZSolveMatrix::TotalStrain(TPZFMatrix<REAL> &delta_strain, TPZFMatrix<REAL> &total_strain) {
    total_strain = total_strain + delta_strain;
}

void TPZSolveMatrix::ElasticStrain(TPZFMatrix<REAL> &total_strain, TPZFMatrix<REAL> &plastic_strain, TPZFMatrix<REAL> &elastic_strain) {
    elastic_strain = total_strain - plastic_strain;
}

void TPZSolveMatrix::PlasticStrain(TPZFMatrix<REAL> &total_strain, TPZFMatrix<REAL> &elastic_strain, TPZFMatrix<REAL> &plastic_strain) {
    plastic_strain = total_strain - elastic_strain;
}

//Compute stress
void TPZSolveMatrix::ComputeStress(TPZFMatrix<REAL> &elastic_strain, TPZFMatrix<REAL> &sigma) {
    REAL lambda = fMaterialData.ElasticResponse().Lambda();
    REAL mu = fMaterialData.ElasticResponse().Mu();
    sigma.Resize(4*fNpts,1);

    for (int64_t ipts=0; ipts < fNpts; ipts++) {
        //plane strain
        sigma(4 * ipts, 0) = elastic_strain(2 * ipts, 0) * (lambda + 2. * mu) + elastic_strain(2 * ipts + 2 * fNpts + 1, 0) * lambda; // Sigma xx
        sigma(4 * ipts + 1, 0) = elastic_strain(2 * ipts + 2 * fNpts + 1, 0) * (lambda + 2. * mu) + elastic_strain(2 * ipts, 0) * lambda; // Sigma yy
        sigma(4 * ipts + 2, 0) = lambda * (elastic_strain(2 * ipts, 0) + elastic_strain(2 * ipts + 2 * fNpts + 1, 0)); // Sigma zz
        sigma(4 * ipts + 3, 0) = mu * (elastic_strain(2 * ipts + 1, 0) + elastic_strain(2 * ipts + 2 * fNpts, 0)); // Sigma xy
    }
}

//Compute strain
void TPZSolveMatrix::ComputeStrain(TPZFMatrix<REAL> &sigma, TPZFMatrix<REAL> &elastic_strain) {
    REAL E = fMaterialData.ElasticResponse().E();
    REAL nu = fMaterialData.ElasticResponse().Poisson();
    int dim = fCmesh->Dimension();

    for (int ipts = 0; ipts < fNpts; ipts++) {
        elastic_strain(2 * ipts + 0, 0) = 1/fWeight[ipts]*(1. / E * (sigma(2*ipts,0)*(1.-nu*nu) - sigma(2*ipts + 2*fNpts +1,0)*(nu+nu*nu))); //exx
        elastic_strain(2 * ipts + 1, 0) = 1/fWeight[ipts]*((1. + nu) / E * sigma(2 * ipts + 1, 0)); //exy
        elastic_strain(2 * ipts + dim * fNpts + 0, 0) = elastic_strain(2 * ipts + 1, 0); //exy
        elastic_strain(2 * ipts + dim * fNpts + 1, 0) = 1/fWeight[ipts]*(1. / E * (sigma(2*ipts + 2*fNpts +1,0)*(1.-nu*nu) - sigma(2*ipts,0)*(nu+nu*nu))); //eyy
    }
}

void TPZSolveMatrix::SpectralDecomposition(TPZFMatrix<REAL> &sigma_trial, TPZFMatrix<REAL> &eigenvalues, TPZFMatrix<REAL> &eigenvectors) {
    REAL maxel;
    TPZVec<REAL> interval(2);
    eigenvalues.Resize(3*fNpts,1);
    eigenvectors.Resize(9*fNpts,1);

    for (int64_t ipts = 0; ipts < fNpts; ipts++) {
        Normalize(&sigma_trial(4*ipts, 0), maxel);
        Interval(&sigma_trial(4*ipts, 0), &interval[0]);
        NewtonIterations(&interval[0], &sigma_trial(4*ipts, 0), &eigenvalues(3*ipts, 0), maxel);
        Eigenvectors(&sigma_trial(4*ipts, 0), &eigenvalues(3*ipts, 0), &eigenvectors(9*ipts,0),maxel);
    }
}

void TPZSolveMatrix::ProjectSigma(TPZFMatrix<REAL> &eigenvalues, TPZFMatrix<REAL> &sigma_projected) {
    int dim = fCmesh->Dimension();

    REAL mc_psi = fMaterialData.FrictionAngle();

    sigma_projected.Resize(3*fNpts,1);
    sigma_projected.Zero();
    TPZFMatrix<REAL> elastic_strain_np1(dim*dim*fNpts);

    TPZFMatrix<REAL> m_type(fNpts, 1, 0.);
    TPZFMatrix<REAL> alpha(fNpts, 1, 0.);
    bool check = false;

    for (int ipts = 0; ipts < fNpts; ipts++) {
        m_type(ipts,0) = 0;
        check = PhiPlane(&eigenvalues(3*ipts, 0), &sigma_projected(3*ipts, 0)); //elastic domain
        if (!check) { //plastic domain
            m_type(ipts,0) = 1;
            check = ReturnMappingMainPlane(&eigenvalues(3*ipts, 0), &sigma_projected(3*ipts, 0), alpha(ipts,0)); //main plane
            if (!check) { //edges or apex
                if  (((1 - sin(mc_psi)) * eigenvalues(0 + 3*ipts, 0) - 2. * eigenvalues(1 + 3*ipts, 0) + (1 + sin(mc_psi)) * eigenvalues(2 + 3*ipts, 0)) > 0) { // right edge
                    check = ReturnMappingRightEdge(&eigenvalues(3*ipts, 0), &sigma_projected(3*ipts, 0), alpha(ipts,0));
                } else { //left edge
                    check = ReturnMappingLeftEdge(&eigenvalues(3*ipts, 0), &sigma_projected(3*ipts, 0), alpha(ipts,0));
                }
                if (!check) { //apex
                    m_type(ipts,0) = -1;
                    ReturnMappingApex(&eigenvalues(3*ipts, 0), &sigma_projected(3*ipts, 0), alpha(ipts,0));
                }
            }
        }
    }
}

void TPZSolveMatrix::StressCompleteTensor(TPZFMatrix<REAL> &sigma_projected, TPZFMatrix<REAL> &eigenvectors, TPZFMatrix<REAL> &sigma){
    sigma.Resize(4*fNpts,1);
    int dim = fCmesh->Dimension();

    for (int ipts = 0; ipts < fNpts; ipts++) {
        sigma(2*ipts + 0,0) = fWeight[ipts]*(sigma_projected(3*ipts + 0,0)*eigenvectors(9*ipts + 0,0)*eigenvectors(9*ipts + 0,0) + sigma_projected(3*ipts + 1,0)*eigenvectors(9*ipts + 3,0)*eigenvectors(9*ipts + 3,0) + sigma_projected(3*ipts + 2,0)*eigenvectors(9*ipts + 6,0)*eigenvectors(9*ipts + 6,0));
        sigma(2*ipts + 1,0) = fWeight[ipts]*(sigma_projected(3*ipts + 0,0)*eigenvectors(9*ipts + 0,0)*eigenvectors(9*ipts + 1,0) + sigma_projected(3*ipts + 1,0)*eigenvectors(9*ipts + 3,0)*eigenvectors(9*ipts + 4,0) + sigma_projected(3*ipts + 2,0)*eigenvectors(9*ipts + 6,0)*eigenvectors(9*ipts + 7,0));
        sigma(2*ipts + dim*fNpts,0) = sigma(2*ipts + 1,0);
        sigma(2*ipts + dim*fNpts + 1,0) = fWeight[ipts]*(sigma_projected(3*ipts + 0,0)*eigenvectors(9*ipts + 1,0)*eigenvectors(9*ipts + 1,0) + sigma_projected(3*ipts + 1,0)*eigenvectors(9*ipts + 4,0)*eigenvectors(9*ipts + 4,0) + sigma_projected(3*ipts + 2,0)*eigenvectors(9*ipts + 7,0)*eigenvectors(9*ipts + 7,0));
    }
}

void TPZSolveMatrix::NodalForces(TPZFMatrix<REAL> &sigma, TPZFMatrix<REAL> &nodal_forces) {
    int64_t nelem = fRowSizes.size();
    int64_t npts = fNpts;
    int dim = fCmesh->Dimension();
    int64_t size = dim*fNpts;
    nodal_forces.Resize(size,1);
    nodal_forces.Zero();

    for (int iel = 0; iel < nelem; iel++) {
        for (int i = 0; i < fColSizes[iel]; i++) {
            for (int k = 0; k < fRowSizes[iel]; k++) {
                nodal_forces(i + fColFirstIndex[iel], 0) += fStorage[k + i * fRowSizes[iel] + fMatrixPosition[iel]] * sigma(k + fRowFirstIndex[iel], 0);
                nodal_forces(i + fColFirstIndex[iel] + size/dim, 0) +=  fStorage[k + i * fRowSizes[iel] + fMatrixPosition[iel]] * sigma(k + fRowFirstIndex[iel] + size, 0);
            }
        }
    }
}

void TPZSolveMatrix::ColoredAssemble(TPZFMatrix<STATE>  &nodal_forces_vec, TPZFMatrix<STATE> &nodal_forces_global) {
    int64_t ncolor = *std::max_element(fElemColor.begin(), fElemColor.end())+1;
    int64_t sz = fIndexes.size();
    int64_t neq = fCmesh->NEquations();
    nodal_forces_global.Resize(neq*ncolor,1);
    nodal_forces_global.Zero();


    cblas_dsctr(sz, nodal_forces_vec, &fIndexesColor[0], &nodal_forces_global(0,0));

    int64_t colorassemb = ncolor / 2.;
    while (colorassemb > 0) {

        int64_t firsteq = (ncolor - colorassemb) * neq;
        cblas_daxpy(firsteq, 1., &nodal_forces_global(firsteq, 0), 1., &nodal_forces_global(0, 0), 1.);

        ncolor -= colorassemb;
        colorassemb = ncolor/2;
    }
    nodal_forces_global.Resize(neq, 1);
}

void TPZSolveMatrix::ColoringElements() const {
    int64_t nelem_c = fCmesh->NElements();
    int64_t nconnects = fCmesh->NConnects();
    TPZVec<int64_t> connects_vec(nconnects,0);

    int64_t contcolor = 0;
    bool needstocontinue = true;

    while (needstocontinue)
    {
        int it = 0;
        needstocontinue = false;
        for (int64_t iel = 0; iel < nelem_c; iel++) {
            TPZCompEl *cel = fCmesh->Element(iel);
            if (!cel || cel->Dimension() != fCmesh->Dimension()) continue;

            it++;
            if (fElemColor[it-1] != -1) continue;

            TPZStack<int64_t> connectlist;
            fCmesh->Element(iel)->BuildConnectList(connectlist);
            int64_t ncon = connectlist.size();

            int64_t icon;
            for (icon = 0; icon < ncon; icon++) {
                if (connects_vec[connectlist[icon]] != 0) break;
            }
            if (icon != ncon) {
                needstocontinue = true;
                continue;
            }
            fElemColor[it-1] = contcolor;
//            cel->Reference()->SetMaterialId(contcolor);

            for (icon = 0; icon < ncon; icon++) {
                connects_vec[connectlist[icon]] = 1;
            }
        }
        contcolor++;
        connects_vec.Fill(0);
    }


    int64_t nelem = fRowSizes.size();
    int64_t neq = fCmesh->NEquations();
    for (int64_t iel = 0; iel < nelem; iel++) {
        int64_t cols = fColSizes[iel];
        int64_t cont_cols = fColFirstIndex[iel];

        for (int64_t icols = 0; icols < cols; icols++) {
            fIndexesColor[cont_cols + icols] = fIndexes[cont_cols + icols] + fElemColor[iel]*neq;
            fIndexesColor[cont_cols+ fNphis + icols] = fIndexes[cont_cols + fNphis + icols] + fElemColor[iel]*neq;
        }
    }
}

TPZFMatrix<REAL> TPZSolveMatrix::AssembleResidual() {
    TPZFMatrix<REAL> gather_solution;
    TPZFMatrix<REAL> delta_strain;
    TPZFMatrix<REAL> elastic_strain;
    TPZFMatrix<REAL> sigma_trial;
    TPZFMatrix<REAL> eigenvalues;
    TPZFMatrix<REAL> eigenvectors;
    TPZFMatrix<REAL> sigma_projected;
    TPZFMatrix<REAL> sigma;
    TPZFMatrix<REAL> nodal_forces;
    TPZFMatrix<REAL> residual;

    //residual assemble
    GatherSolution(fSolution, gather_solution);
    DeltaStrain(gather_solution, delta_strain);
    TotalStrain(delta_strain, fTotalStrain);
    ElasticStrain(fTotalStrain, fPlasticStrain, elastic_strain);
    ComputeStress(elastic_strain, sigma_trial);
    SpectralDecomposition(sigma_trial, eigenvalues, eigenvectors);
    ProjectSigma(eigenvalues, sigma_projected);
    StressCompleteTensor(sigma_projected, eigenvectors, sigma);
    NodalForces(sigma, nodal_forces);
    ColoredAssemble(nodal_forces,residual);

    //update strain
    ComputeStrain(sigma, elastic_strain);
    PlasticStrain(fTotalStrain, elastic_strain, fPlasticStrain);

    return residual;
}

void TPZSolveMatrix::SetDataStructure(){
    int dim_mesh = (fCmesh->Reference())->Dimension(); // Mesh dimension
    int64_t nelem_c = fCmesh->NElements(); // Number of computational elements
    std::vector<int64_t> cel_indexes;

// Number of domain geometric elements
    for (int64_t i = 0; i < nelem_c; i++) {
        TPZCompEl *cel = fCmesh->Element(i);
        if (!cel) continue;
        TPZGeoEl *gel = fCmesh->Element(i)->Reference();
        if (!gel || gel->Dimension() != dim_mesh) continue;
        cel_indexes.push_back(cel->Index());
    }

    if (cel_indexes.size() == 0) {
        DebugStop();
    }

// RowSizes and ColSizes vectors
    int64_t nelem = cel_indexes.size();
    TPZVec<MKL_INT> rowsizes(nelem);
    TPZVec<MKL_INT> colsizes(nelem);

    int64_t npts_tot = 0;
    int64_t nf_tot = 0;
    int it = 0;
    for (auto iel : cel_indexes) {
        //Verification
        TPZCompEl *cel = fCmesh->Element(iel);

        //Integration rule
        TPZInterpolatedElement *cel_inter = dynamic_cast<TPZInterpolatedElement * >(cel);
        if (!cel_inter) DebugStop();
        TPZIntPoints *int_rule = &(cel_inter->GetIntegrationRule());

        int64_t npts = int_rule->NPoints(); // number of integration points of the element
        int64_t dim = cel_inter->Dimension(); //dimension of the element
        int64_t nf = cel_inter->NShapeF(); // number of shape functions of the element

        rowsizes[it] = dim * npts;
        colsizes[it] = nf;

        it++;

        npts_tot += npts;
        nf_tot += nf;
    }
    this->SetNumberofIntPoints(npts_tot);
    this->SetNumberofPhis(nf_tot);
    this->SetRowandColSizes(rowsizes, colsizes);

// Dphi matrix, weight and indexes vectors
    TPZFMatrix<REAL> elmatrix;
    TPZStack<REAL> weight;
    TPZManVector<MKL_INT> indexes(dim_mesh * nf_tot);

    int64_t cont1 = 0;
    int64_t cont2 = 0;
    it = 0;
    for (auto iel : cel_indexes) {
        //Verification
        TPZCompEl *cel = fCmesh->Element(iel);

        //Integration rule
        TPZInterpolatedElement *cel_inter = dynamic_cast<TPZInterpolatedElement * >(cel);
        if (!cel_inter) DebugStop();
        TPZIntPoints *int_rule = &(cel_inter->GetIntegrationRule());

        int64_t npts = int_rule->NPoints(); // number of integration points of the element
        int64_t dim = cel_inter->Dimension(); //dimension of the element
        int64_t nf = cel_inter->NShapeF(); // number of shape functions of the element

        TPZMaterialData data;
        cel_inter->InitMaterialData(data);

        elmatrix.Resize(dim * npts, nf);
        for (int64_t inpts = 0; inpts < npts; inpts++) {
            TPZManVector<REAL> qsi(dim, 1);
            REAL w;
            int_rule->Point(inpts, qsi, w);
            cel_inter->ComputeRequiredData(data, qsi);
            weight.Push(w * std::abs(data.detjac)); //weight = w * detjac

            TPZFMatrix<REAL> axes = data.axes;
            TPZFMatrix<REAL> dphix = data.dphix;
            TPZFMatrix<REAL> dphiXY;
            axes.Transpose();
            axes.Multiply(dphix,dphiXY);

            for (int inf = 0; inf < nf; inf++) {
                for (int idim = 0; idim < dim; idim++)
                    elmatrix(inpts * dim + idim, inf) = dphiXY(idim, inf);
            }
        }
        this->SetElementMatrix(it, elmatrix);
        it++;

        //Indexes vector
        int64_t ncon = cel->NConnects();
        for (int64_t icon = 0; icon < ncon; icon++) {
            int64_t id = cel->ConnectIndex(icon);
            TPZConnect &df = fCmesh->ConnectVec()[id];
            int64_t conid = df.SequenceNumber();
            if (df.NElConnected() == 0 || conid < 0 || fCmesh->Block().Size(conid) == 0) continue;
            else {
                int64_t pos = fCmesh->Block().Position(conid);
                int64_t nsize = fCmesh->Block().Size(conid);
                for (int64_t isize = 0; isize < nsize; isize++) {
                    if (isize % 2 == 0) {
                        indexes[cont1] = pos + isize;
                        cont1++;
                    } else {
                        indexes[cont2 + nf_tot] = pos + isize;
                        cont2++;
                    }
                }
            }
        }
    }
    this->SetIndexes(indexes);
    this->SetWeightVector(weight);
    this->ColoringElements();
}
