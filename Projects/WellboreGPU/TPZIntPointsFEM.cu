#include "TPZIntPointsFEM.h"
#include "TPZTensor.h"
#include "pzmatrix.h"
#include <stdlib.h>
#include "TPZTensor.h"
#include "TPZVTKGeoMesh.h"
#include "pzintel.h"
#include "pzskylstrmatrix.h"

#define NT 512

__device__ void Normalize(REAL *sigma, REAL &maxel) {
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

__device__ void Interval(REAL *sigma, REAL *interval) {
	__shared__ REAL lower_vec[3];
	__shared__ REAL upper_vec[3];

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

__device__ void NewtonIterations(REAL *interval, REAL *sigma, REAL *eigenvalues, REAL &maxel) {
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

	//sorting in descending order
	for (int i = 0; i < 3; ++i) {
		for (int j = i + 1; j < 3; ++j) {
			if (eigenvalues[i] < eigenvalues[j]) {
				REAL a = eigenvalues[i];
				eigenvalues[i] = eigenvalues[j];
				eigenvalues[j] = a;
			}
		}
	}
}

__device__ void Multiplicity1(REAL *sigma, REAL eigenvalue, REAL *eigenvector) {
	__shared__ REAL det[3];
	det[0] = (sigma[0] - eigenvalue) * (sigma[1] - eigenvalue) - sigma[3] * sigma[3];
	det[1] = (sigma[0] - eigenvalue) * (sigma[2] - eigenvalue);
	det[2] = (sigma[1] - eigenvalue) * (sigma[2] - eigenvalue);

	REAL maxdet = fabs(det[0]);
	for (int i = 1; i < 3; i++) {
		if (fabs(det[i]) > fabs(maxdet)) {
			maxdet = fabs(det[i]);
		}
	}
	__shared__ REAL v[3];
	if (maxdet == fabs(det[0])) {
		v[0] = 0;
		v[1] = 0;
		v[2] = 1;

	} else if (maxdet == fabs(det[1])) {
		v[0] = 1 / det[1] * (-(sigma[2] - eigenvalue) * sigma[3]);
		v[1] = 1;
		v[2] = 0;

	} else {
		v[0] = 1;
		v[1] = 1 / det[2] * (-(sigma[2] - eigenvalue) * sigma[3]);
		v[2] = 0;
	}
	REAL norm = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
	eigenvector[0] = v[0] / norm;
	eigenvector[1] = v[1] / norm;
	eigenvector[2] = v[2] / norm;
}

__device__ void Multiplicity2(REAL *sigma, REAL eigenvalue, REAL *eigenvector1,
		REAL *eigenvector2) {
	__shared__ REAL x[3];
	x[0] = sigma[0] - eigenvalue;
	x[1] = sigma[1] - eigenvalue;
	x[2] = sigma[2] - eigenvalue;

	REAL maxx = fabs(x[0]);
	for (int i = 1; i < 3; i++) {
		if (fabs(x[i]) > fabs(maxx)) {
			maxx = fabs(x[i]);
		}
	}

	__shared__ REAL v1[3];
	__shared__ REAL v2[3];

	if (maxx == fabs(x[0])) {
		v1[0] = -sigma[3] / x[0];
		v1[1] = 1;
		v1[2] = 0;

		v2[0] = 0;
		v2[1] = 0;
		v2[2] = 1;

	} else if (maxx == fabs(x[1])) {
		v1[0] = 1;
		v1[1] = -sigma[3] / x[1];
		v1[2] = 0;

		v2[0] = 0;
		v2[1] = 0;
		v2[2] = 1;

	} else {
		v1[0] = 1;
		v1[1] = 0;
		v1[2] = 0;

		v2[0] = 0;
		v2[1] = 1;
		v2[2] = 0;

	}
	REAL norm1 = sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]);
	REAL norm2 = sqrt(v2[0] * v2[0] + v2[1] * v1[1] + v2[2] * v2[2]);

	eigenvector1[0] = v1[0] / norm1;
	eigenvector1[1] = v1[1] / norm1;
	eigenvector1[2] = v1[2] / norm1;

	eigenvector2[0] = v2[0] / norm2;
	eigenvector2[1] = v2[1] / norm2;
	eigenvector2[2] = v2[2] / norm2;
}

__device__ void Eigenvectors(REAL *sigma, REAL *eigenvalues, REAL *eigenvectors,
		REAL &maxel) {
	sigma[0] *= maxel;
	sigma[1] *= maxel;
	sigma[2] *= maxel;
	sigma[3] *= maxel;

	if ((eigenvalues[0] == eigenvalues[1])
			&& (eigenvalues[1] == eigenvalues[2])) {
		eigenvectors[0] = 1.;
		eigenvectors[1] = 0.;
		eigenvectors[2] = 0.;

		eigenvectors[3] = 0.;
		eigenvectors[4] = 1.;
		eigenvectors[5] = 0.;

		eigenvectors[6] = 0.;
		eigenvectors[7] = 0.;
		eigenvectors[8] = 1.;
	} else {
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

__global__ void DeltaStrainKernel(int64_t nelem, REAL *storage, int *rowsizes, int *colsizes, int *matrixpos, int *rowfirstindex, int* colfirstindex, int npts, int nphis, REAL *gather_solution, REAL *delta_strain) {
	int iel = blockIdx.x * blockDim.x + threadIdx.x;

	if (iel < nelem) {
		for (int i = 0; i < rowsizes[iel]; i++) {
			for (int k = 0; k < colsizes[iel]; k++) {
				delta_strain[i + rowfirstindex[iel]] += storage[k * rowsizes[iel] + i + matrixpos[iel]] * gather_solution[k + colfirstindex[iel]];
				delta_strain[i + rowfirstindex[iel] + npts] += storage[k * rowsizes[iel] + i + matrixpos[iel]] * gather_solution[k + colfirstindex[iel] + nphis];
			}
		}
	}
}

__global__ void NodalForcesKernel(int64_t nelem, REAL *storage, int *rowsizes, int *colsizes, int *matrixpos, int *rowfirstindex, int* colfirstindex, int npts, int nphis, REAL *sigma, REAL *nodal_forces) {
	int iel = blockIdx.x * blockDim.x + threadIdx.x;

	if (iel < nelem) {
		for (int i = 0; i < colsizes[iel]; i++) {
			for (int k = 0; k < rowsizes[iel]; k++) {
				nodal_forces[i + colfirstindex[iel]] -= storage[k + i * rowsizes[iel] + matrixpos[iel]] * sigma[k + rowfirstindex[iel]];
				nodal_forces[i + colfirstindex[iel] + nphis] -= storage[k + i * rowsizes[iel] + matrixpos[iel]] * sigma[k + rowfirstindex[iel] + npts];
			}
		}
	}
}

__global__ void ComputeStressKernel(int64_t fNpts, int fDim,
		REAL *elastic_strain, REAL *sigma, REAL mu, REAL lambda) {
	int ipts = blockIdx.x * blockDim.x + threadIdx.x;

	if (ipts < fNpts / fDim) {
		sigma[4 * ipts] = elastic_strain[2 * ipts] * (lambda + 2. * mu) + elastic_strain[2 * ipts + fNpts + 1] * lambda; // Sigma xx
		sigma[4 * ipts + 1] = elastic_strain[2 * ipts + fNpts + 1] * (lambda + 2. * mu) + elastic_strain[2 * ipts] * lambda; // Sigma yy
		sigma[4 * ipts + 2] = lambda * (elastic_strain[2 * ipts] + elastic_strain[2 * ipts + fNpts + 1]); // Sigma zz
		sigma[4 * ipts + 3] = mu * (elastic_strain[2 * ipts + 1] + elastic_strain[2 * ipts + fNpts]); // Sigma xy
	}
}

__global__ void ComputeStrainKernel(int64_t fNpts, int fDim, REAL *sigma, REAL *elastic_strain, REAL nu, REAL E, REAL *weight) {
	int ipts = blockIdx.x * blockDim.x + threadIdx.x;

	if (ipts < fNpts / fDim) {
		elastic_strain[2 * ipts + 0] = 1 / weight[ipts] * (1. / E * (sigma[2 * ipts] * (1. - nu * nu) - sigma[2 * ipts + fNpts + 1] * (nu + nu * nu))); //exx
		elastic_strain[2 * ipts + 1] = 1 / weight[ipts] * ((1. + nu) / E * sigma[2 * ipts + 1]); //exy
		elastic_strain[2 * ipts + fNpts + 0] = elastic_strain[2 * ipts + 1]; //exy
		elastic_strain[2 * ipts + fNpts + 1] = 1 / weight[ipts] * (1. / E * (sigma[2 * ipts + fNpts + 1] * (1. - nu * nu) - sigma[2 * ipts] * (nu + nu * nu))); //eyy
	}
}

__global__ void SpectralDecompositionKernel(int64_t fNpts, int fDim, REAL *sigma_trial, REAL *eigenvalues, REAL *eigenvectors) {
	int ipts = blockIdx.x * blockDim.x + threadIdx.x;

	__shared__ REAL maxel;
	__shared__ REAL interval[2];
	if (ipts < fNpts / fDim) {
		Normalize(&sigma_trial[4 * ipts], maxel);
		Interval(&sigma_trial[4 * ipts], &interval[0]);
		NewtonIterations(&interval[0], &sigma_trial[4 * ipts], &eigenvalues[3 * ipts], maxel);
		Eigenvectors(&sigma_trial[4 * ipts], &eigenvalues[3 * ipts], &eigenvectors[9 * ipts], maxel);
	}
}

__device__ bool PhiPlane(REAL *eigenvalues, REAL *sigma_projected, REAL mc_phi, REAL mc_cohesion) {
	const REAL sinphi = sin(mc_phi);
	const REAL cosphi = cos(mc_phi);

	REAL phi = eigenvalues[0] - eigenvalues[2]
			+ (eigenvalues[0] + eigenvalues[2]) * sinphi
			- 2. * mc_cohesion * cosphi;

	sigma_projected[0] = eigenvalues[0];
	sigma_projected[1] = eigenvalues[1];
	sigma_projected[2] = eigenvalues[2];

	bool check_validity = (fabs(phi) < 1.e-12) || (phi < 0.0);
	return check_validity;
}

__device__ bool ReturnMappingMainPlane(REAL *eigenvalues, REAL *sigma_projected, REAL &m_hardening, REAL mc_phi, REAL mc_psi, REAL mc_cohesion, REAL K, REAL G) {
	const REAL sinphi = sin(mc_phi);
	const REAL sinpsi = sin(mc_psi);
	const REAL cosphi = cos(mc_phi);
	const REAL sinphi2 = sinphi * sinphi;
	const REAL cosphi2 = 1. - sinphi2;
	const REAL constA = 4. * G * (1. + sinphi * sinpsi / 3.) + 4. * K * sinphi * sinpsi;

	REAL phi = eigenvalues[0] - eigenvalues[2] + (eigenvalues[0] + eigenvalues[2]) * sinphi - 2. * mc_cohesion * cosphi;

	REAL gamma = 0;
	int n_iterations = 30;
	for (int i = 0; i < n_iterations; i++) {
		double jac = -constA - 4. * cosphi2 * 0; // H=0
		double delta_gamma = -phi / jac;
		gamma += delta_gamma;
		phi = eigenvalues[0] - eigenvalues[2] + (eigenvalues[0] + eigenvalues[2]) * sinphi - 2. * mc_cohesion * cosphi - constA * gamma;
		if (fabs(phi) < 1.e-12) {
			break;
		}
	}

	eigenvalues[0] -= (2. * G * (1 + sinpsi / 3.) + 2. * K * sinpsi) * gamma;
	eigenvalues[1] += (4. * G / 3. - K * 2.) * sinpsi * gamma;
	eigenvalues[2] += (2. * G * (1 - sinpsi / 3.) - 2. * K * sinpsi) * gamma;
	sigma_projected[0] = eigenvalues[0];
	sigma_projected[1] = eigenvalues[1];
	sigma_projected[2] = eigenvalues[2];

	m_hardening += gamma * 2. * cosphi;

	bool check_validity = (eigenvalues[0] > eigenvalues[1] || fabs(eigenvalues[0] - eigenvalues[1]) < 1.e-12) && (eigenvalues[1] > eigenvalues[2] || fabs(eigenvalues[1] - eigenvalues[2]) < 1.e-12);
	return check_validity;
}

__device__ bool ReturnMappingRightEdge(REAL *eigenvalues, REAL *sigma_projected, REAL &m_hardening, REAL mc_phi, REAL mc_psi, REAL mc_cohesion, REAL K, REAL G) {
	const REAL sinphi = sin(mc_phi);
	const REAL sinpsi = sin(mc_psi);
	const REAL cosphi = cos(mc_phi);

	__shared__ REAL gamma[2], phi[2], sigma_bar[2], ab[2];

	__shared__ REAL jac[2][2], jac_inv[2][2];

	sigma_bar[0] = eigenvalues[0] - eigenvalues[2] + (eigenvalues[0] + eigenvalues[2]) * sinphi;
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

		phi[0] = sigma_bar[0] - ab[0] * gamma[0] - ab[1] * gamma[1]
				- 2. * cosphi * mc_cohesion;
		phi[1] = sigma_bar[1] - ab[1] * gamma[0] - ab[0] * gamma[1]
				- 2. * cosphi * mc_cohesion;

		double res = (fabs(phi[0]) + fabs(phi[1]));

		if (fabs(res) < 1.e-12) {
			break;
		}
	}

	eigenvalues[0] -= (2. * G * (1 + sinpsi / 3.) + 2. * K * sinpsi) * (gamma[0] + gamma[1]);
	eigenvalues[1] += ((4. * G / 3. - K * 2.) * sinpsi) * gamma[0] + (2. * G * (1. - sinpsi / 3.) - 2. * K * sinpsi) * gamma[1];
	eigenvalues[2] += (2. * G * (1 - sinpsi / 3.) - 2. * K * sinpsi) * gamma[0] + ((4. * G / 3. - 2. * K) * sinpsi) * gamma[1];
	sigma_projected[0] = eigenvalues[0];
	sigma_projected[1] = eigenvalues[1];
	sigma_projected[2] = eigenvalues[2];

	m_hardening += (gamma[0] + gamma[1]) * 2. * cosphi;

	bool check_validity = (eigenvalues[0] > eigenvalues[1] || fabs(eigenvalues[0] - eigenvalues[1]) < 1.e-12)
			&& (eigenvalues[1] > eigenvalues[2] || fabs(eigenvalues[1] - eigenvalues[2]) < 1.e-12);
	return check_validity;
}

__device__ bool ReturnMappingLeftEdge(REAL *eigenvalues, REAL *sigma_projected, REAL &m_hardening, REAL mc_phi, REAL mc_psi, REAL mc_cohesion, REAL K, REAL G) {
	const REAL sinphi = sin(mc_phi);
	const REAL sinpsi = sin(mc_psi);
	const REAL cosphi = cos(mc_phi);
	const REAL sinphi2 = sinphi * sinphi;
	const REAL cosphi2 = 1. - sinphi2;

	__shared__ REAL gamma[2], phi[2], sigma_bar[2], ab[2];

	__shared__ REAL jac[2][2], jac_inv[2][2];

	sigma_bar[0] = eigenvalues[0] - eigenvalues[2] + (eigenvalues[0] + eigenvalues[2]) * sinphi;
	sigma_bar[1] = eigenvalues[1] - eigenvalues[2] + (eigenvalues[1] + eigenvalues[2]) * sinphi;

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

	eigenvalues[0] += -(2. * G * (1 + sinpsi / 3.) + 2. * K * sinpsi) * gamma[0] + ((4. * G / 3. - 2. * K) * sinpsi) * gamma[1];
	eigenvalues[1] += ((4. * G / 3. - K * 2.) * sinpsi) * gamma[0] - (2. * G * (1. + sinpsi / 3.) + 2. * K * sinpsi) * gamma[1];
	eigenvalues[2] += (2. * G * (1 - sinpsi / 3.) - 2. * K * sinpsi) * (gamma[0] + gamma[1]);
	sigma_projected[0] = eigenvalues[0];
	sigma_projected[1] = eigenvalues[1];
	sigma_projected[2] = eigenvalues[2];

	m_hardening += (gamma[0] + gamma[1]) * 2. * cosphi;

	bool check_validity = (eigenvalues[0] > eigenvalues[1] || fabs(eigenvalues[0] - eigenvalues[1]) < 1.e-12) && (eigenvalues[1] > eigenvalues[2] || fabs(eigenvalues[1] - eigenvalues[2]) < 1.e-12);
	return check_validity;
}

__device__ bool ReturnMappingApex(REAL *eigenvalues, REAL *sigma_projected, REAL &m_hardening, REAL mc_phi, REAL mc_psi, REAL mc_cohesion, REAL K) {
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

__global__ void ProjectSigmaKernel(int64_t fNpts, int fDim, REAL mc_phi, REAL mc_psi, REAL mc_cohesion, REAL K, REAL G, REAL *eigenvalues, REAL *sigma_projected, REAL *m_type, REAL *alpha) {
	int ipts = blockIdx.x * blockDim.x + threadIdx.x;

	bool check = false;
	if (ipts < fNpts / fDim) {
		m_type[ipts] = 0;
		check = PhiPlane(&eigenvalues[3 * ipts], &sigma_projected[3 * ipts], mc_phi, mc_cohesion); //elastic domain
		if (!check) { //plastic domain
			m_type[ipts] = 1;
			check = ReturnMappingMainPlane(&eigenvalues[3 * ipts], &sigma_projected[3 * ipts], alpha[ipts], mc_phi, mc_psi, mc_cohesion, K, G); //main plane
			if (!check) { //edges or apex
				if (((1 - sin(mc_psi)) * eigenvalues[0 + 3 * ipts] - 2. * eigenvalues[1 + 3 * ipts] + (1 + sin(mc_psi)) * eigenvalues[2 + 3 * ipts]) > 0) { // right edge
					check = ReturnMappingRightEdge(&eigenvalues[3 * ipts], &sigma_projected[3 * ipts], alpha[ipts], mc_phi, mc_psi, mc_cohesion, K, G);
				} else { //left edge
					check = ReturnMappingLeftEdge(&eigenvalues[3 * ipts], &sigma_projected[3 * ipts], alpha[ipts], mc_phi, mc_psi, mc_cohesion, K, G);
				}
				if (!check) { //apex
					m_type[ipts] = -1;
					ReturnMappingApex(&eigenvalues[3 * ipts], &sigma_projected[3 * ipts], alpha[ipts], mc_phi, mc_psi, mc_cohesion, K);
				}
			}
		}
	}
}

__global__ void StressCompleteTensorKernel(int64_t fNpts, int fDim, REAL *sigma_projected, REAL *eigenvectors, REAL *sigma, REAL *weight) {
	int ipts = blockIdx.x * blockDim.x + threadIdx.x;

	if (ipts < fNpts / fDim) {

		sigma[2 * ipts + 0] = weight[ipts] * (sigma_projected[3 * ipts + 0] * eigenvectors[9 * ipts + 0] * eigenvectors[9 * ipts + 0] + sigma_projected[3 * ipts + 1]
								* eigenvectors[9 * ipts + 3] * eigenvectors[9 * ipts + 3] + sigma_projected[3 * ipts + 2] * eigenvectors[9 * ipts + 6] * eigenvectors[9 * ipts + 6]);
		sigma[2 * ipts + 1] = weight[ipts] * (sigma_projected[3 * ipts + 0] * eigenvectors[9 * ipts + 0] * eigenvectors[9 * ipts + 1] + sigma_projected[3 * ipts + 1]
								* eigenvectors[9 * ipts + 3] * eigenvectors[9 * ipts + 4] + sigma_projected[3 * ipts + 2] * eigenvectors[9 * ipts + 6] * eigenvectors[9 * ipts + 7]);
		sigma[2 * ipts + fNpts] = sigma[2 * ipts + 1];
		sigma[2 * ipts + fNpts + 1] = weight[ipts] * (sigma_projected[3 * ipts + 0] * eigenvectors[9 * ipts + 1] * eigenvectors[9 * ipts + 1] + sigma_projected[3 * ipts + 1]
								* eigenvectors[9 * ipts + 4] * eigenvectors[9 * ipts + 4] + sigma_projected[3 * ipts + 2] * eigenvectors[9 * ipts + 7] * eigenvectors[9 * ipts + 7]);
	}
}


TPZIntPointsFEM::TPZIntPointsFEM() :
		fDim(-1), fBoundaryElements(), fCmesh(0), fNpts(-1), fNphis(-1), fElemColor(
				0), fMaterial(0), fRhs(0, 0), fRhsBoundary(0, 0), fSolution(0,
				0), fPlasticStrain(0, 0), fStorage(0), fRowSizes(0), fColSizes(
				0), fMatrixPosition(0), fRowFirstIndex(0), fColFirstIndex(0), fIndexes(
				0), fIndexesColor(0), fWeight() {
//	handle_cusparse = new cusparseHandle_t;
//	handle_cublas = new cublasHandle_t;

	dRhs = new REAL[0];
	dRhsBoundary = new REAL[0];
	dSolution = new REAL[0];
	dPlasticStrain = new REAL[0];
	dStorage = new REAL[0];
	dRowSizes = new int[0];
	dColSizes = new int[0];
	dMatrixPosition = new int[0];
	dRowFirstIndex = new int[0];
	dColFirstIndex = new int[0];
	dIndexes = new int[0];
	dIndexesColor = new int[0];
	dWeight = new REAL[0];

}

TPZIntPointsFEM::TPZIntPointsFEM(TPZCompMesh *cmesh, int materialid) :
		fDim(-1), fBoundaryElements(), fCmesh(0), fNpts(-1), fNphis(-1), fElemColor(
				0), fMaterial(0), fRhs(0, 0), fRhsBoundary(0, 0), fSolution(0,
				0), fPlasticStrain(0, 0), fStorage(0), fRowSizes(0), fColSizes(
				0), fMatrixPosition(0), fRowFirstIndex(0), fColFirstIndex(0), fIndexes(
				0), fIndexesColor(0), fWeight() {
	SetCompMesh(cmesh);
	SetMaterialId(materialid);
//	handle_cusparse = new cusparseHandle_t;
//	handle_cublas = new cublasHandle_t;
	dRhs = new REAL[0];
	dRhsBoundary = new REAL[0];
	dSolution = new REAL[0];
	dPlasticStrain = new REAL[0];
	dStorage = new REAL[0];
	dRowSizes = new int[0];
	dColSizes = new int[0];
	dMatrixPosition = new int[0];
	dRowFirstIndex = new int[0];
	dColFirstIndex = new int[0];
	dIndexes = new int[0];
	dIndexesColor = new int[0];
	dWeight = new REAL[0];
}

TPZIntPointsFEM::~TPZIntPointsFEM() {
	cudaFree(dRhs);
	cudaFree(dRhsBoundary);
	cudaFree(dSolution);
	cudaFree(dPlasticStrain);
	cudaFree(dStorage);
	cudaFree(dRowSizes);
	cudaFree(dColSizes);
	cudaFree(dMatrixPosition);
	cudaFree(dRowFirstIndex);
	cudaFree(dColFirstIndex);
	cudaFree(dIndexes);
	cudaFree(dIndexesColor);
	cudaFree(dWeight);

	cublasDestroy(handle_cublas);
	cusparseDestroy(handle_cusparse);
}

TPZIntPointsFEM::TPZIntPointsFEM(const TPZIntPointsFEM &copy) {
	fDim = copy.fDim;
	fBoundaryElements = copy.fBoundaryElements;
	fCmesh = copy.fCmesh;
	fNpts = copy.fNpts;
	fNphis = copy.fNphis;
	fElemColor = copy.fElemColor;
	fMaterial = copy.fMaterial;

	fRhs = copy.fRhs;
	fRhsBoundary = copy.fRhsBoundary;
	fSolution = copy.fSolution;
	fPlasticStrain = copy.fPlasticStrain;
	fStorage = copy.fStorage;
	fColSizes = copy.fColSizes;
	fRowSizes = copy.fRowSizes;
	fMatrixPosition = copy.fMatrixPosition;
	fRowFirstIndex = copy.fRowFirstIndex;
	fColFirstIndex = copy.fColFirstIndex;
	fIndexes = copy.fIndexes;
	fIndexesColor = copy.fIndexesColor;
	fWeight = copy.fWeight;

	handle_cusparse = copy.handle_cusparse;
	handle_cublas = copy.handle_cublas;

	dRhs = copy.dRhs;
	dRhsBoundary = copy.dRhsBoundary;
	dSolution = copy.dSolution;
	dPlasticStrain = copy.dPlasticStrain;
	dStorage = copy.dStorage;
	dRowSizes = copy.dRowSizes;
	dColSizes = copy.dColSizes;
	dMatrixPosition = copy.dMatrixPosition;
	dRowFirstIndex = copy.dRowFirstIndex;
	dColFirstIndex = copy.dColFirstIndex;
	dIndexes = copy.dIndexes;
	dIndexesColor = copy.dIndexesColor;
	dWeight = copy.dWeight;
}

TPZIntPointsFEM &TPZIntPointsFEM::operator=(const TPZIntPointsFEM &copy) {
	if (&copy == this) {
		return *this;
	}

	fDim = copy.fDim;
	fBoundaryElements = copy.fBoundaryElements;
	fCmesh = copy.fCmesh;
	fNpts = copy.fNpts;
	fNphis = copy.fNphis;
	fElemColor = copy.fElemColor;
	fMaterial = copy.fMaterial;

	fRhs = copy.fRhs;
	fRhsBoundary = copy.fRhsBoundary;
	fSolution = copy.fSolution;
	fPlasticStrain = copy.fPlasticStrain;
	fStorage = copy.fStorage;
	fColSizes = copy.fColSizes;
	fRowSizes = copy.fRowSizes;
	fMatrixPosition = copy.fMatrixPosition;
	fRowFirstIndex = copy.fRowFirstIndex;
	fColFirstIndex = copy.fColFirstIndex;
	fIndexes = copy.fIndexes;
	fIndexesColor = copy.fIndexesColor;
	fWeight = copy.fWeight;

	handle_cusparse = copy.handle_cusparse;
	handle_cublas = copy.handle_cublas;

	dRhs = copy.dRhs;
	dRhsBoundary = copy.dRhsBoundary;
	dSolution = copy.dSolution;
	dPlasticStrain = copy.dPlasticStrain;
	dStorage = copy.dStorage;
	dRowSizes = copy.dRowSizes;
	dColSizes = copy.dColSizes;
	dMatrixPosition = copy.dMatrixPosition;
	dRowFirstIndex = copy.dRowFirstIndex;
	dColFirstIndex = copy.dColFirstIndex;
	dIndexes = copy.dIndexes;
	dIndexesColor = copy.dIndexesColor;
	dWeight = copy.dWeight;

	return *this;
}

void TPZIntPointsFEM::SetDataStructure() {

	int dim_mesh = (fCmesh->Reference())->Dimension(); // Mesh dimension
	this->SetMeshDimension(dim_mesh);
	int64_t nelem_c = fCmesh->NElements(); // Number of computational elements
	std::vector<int64_t> cel_indexes;

// Number of domain geometric elements
	for (int64_t i = 0; i < nelem_c; i++) {
		TPZCompEl *cel = fCmesh->Element(i);
		if (!cel)
			continue;
		TPZGeoEl *gel = fCmesh->Element(i)->Reference();
		if (!gel)
			continue;
		if (gel->Dimension() == dim_mesh)
			cel_indexes.push_back(cel->Index());
		if (gel->Dimension() < dim_mesh)
			fBoundaryElements.Push(cel->Index());
	}

	if (cel_indexes.size() == 0) {
		DebugStop();
	}

// RowSizes and ColSizes vectors
	int64_t nelem = cel_indexes.size();
	TPZVec < MKL_INT > rowsizes(nelem);
	TPZVec < MKL_INT > colsizes(nelem);

	int64_t npts_tot = 0;
	int64_t nf_tot = 0;
	int it = 0;
	for (auto iel : cel_indexes) {
		//Verification
		TPZCompEl *cel = fCmesh->Element(iel);

		//Integration rule
		TPZInterpolatedElement *cel_inter =
				dynamic_cast<TPZInterpolatedElement *>(cel);
		if (!cel_inter)
			DebugStop();
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
	this->SetNumberofIntPoints(dim_mesh * npts_tot);
	this->SetNumberofPhis(nf_tot);
	this->SetRowandColSizes(rowsizes, colsizes);

// Dphi matrix, weight and indexes vectors
	TPZFMatrix < REAL > elmatrix;
	TPZVec < REAL > weight(npts_tot);
	TPZManVector < MKL_INT > indexes(dim_mesh * nf_tot);

	int64_t cont1 = 0;
	int64_t cont2 = 0;
	it = 0;
	int64_t contw = 0;
	for (auto iel : cel_indexes) {
		//Verification
		TPZCompEl *cel = fCmesh->Element(iel);

		//Integration rule
		TPZInterpolatedElement *cel_inter =
				dynamic_cast<TPZInterpolatedElement *>(cel);
		if (!cel_inter)
			DebugStop();
		TPZIntPoints *int_rule = &(cel_inter->GetIntegrationRule());

		int64_t npts = int_rule->NPoints(); // number of integration points of the element
		int64_t dim = cel_inter->Dimension(); //dimension of the element
		int64_t nf = cel_inter->NShapeF(); // number of shape functions of the element

		TPZMaterialData data;
		cel_inter->InitMaterialData(data);

		elmatrix.Resize(dim * npts, nf);
		for (int64_t inpts = 0; inpts < npts; inpts++) {
			TPZManVector < REAL > qsi(dim, 1);
			REAL w;
			int_rule->Point(inpts, qsi, w);
			cel_inter->ComputeRequiredData(data, qsi);
//			weight.Push(w * std::abs(data.detjac)); //weight = w * detjac
			weight[contw] = w * std::abs(data.detjac);
			contw++;

			TPZFMatrix < REAL > axes = data.axes;
			TPZFMatrix < REAL > dphix = data.dphix;
			TPZFMatrix < REAL > dphiXY;
			axes.Transpose();
			axes.Multiply(dphix, dphiXY);

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
			if (df.NElConnected() == 0 || conid < 0
					|| fCmesh->Block().Size(conid) == 0)
				continue;
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
	this->AssembleRhsBoundary();
	this->TransferDataStructure();
}

void TPZIntPointsFEM::ColoringElements() const {
	int64_t nelem_c = fCmesh->NElements();
	int64_t nconnects = fCmesh->NConnects();
	TPZVec<int64_t> connects_vec(nconnects, 0);

	int64_t contcolor = 0;
	bool needstocontinue = true;

	while (needstocontinue) {
		int it = 0;
		needstocontinue = false;
		for (int64_t iel = 0; iel < nelem_c; iel++) {
			TPZCompEl *cel = fCmesh->Element(iel);
			if (!cel || cel->Dimension() != fCmesh->Dimension())
				continue;

			it++;
			if (fElemColor[it - 1] != -1)
				continue;

			TPZStack<int64_t> connectlist;
			fCmesh->Element(iel)->BuildConnectList(connectlist);
			int64_t ncon = connectlist.size();

			int64_t icon;
			for (icon = 0; icon < ncon; icon++) {
				if (connects_vec[connectlist[icon]] != 0)
					break;
			}
			if (icon != ncon) {
				needstocontinue = true;
				continue;
			}
			fElemColor[it - 1] = contcolor;
//            cel->Reference()->SetMaterialId(contcolor);

			for (icon = 0; icon < ncon; icon++) {
				connects_vec[connectlist[icon]] = 1;
			}
		}
		contcolor++;
		connects_vec.Fill(0);
	}
//    ofstream file("colored.vtk");
//    TPZVTKGeoMesh::PrintGMeshVTK(fCmesh->Reference(),file);

	int64_t nelem = fRowSizes.size();
	int64_t neq = fCmesh->NEquations();
	for (int64_t iel = 0; iel < nelem; iel++) {
		int64_t cols = fColSizes[iel];
		int64_t cont_cols = fColFirstIndex[iel];

		for (int64_t icols = 0; icols < cols; icols++) {
			fIndexesColor[cont_cols + icols] = fIndexes[cont_cols + icols]
					+ fElemColor[iel] * neq;
			fIndexesColor[cont_cols + fNphis + icols] = fIndexes[cont_cols
					+ fNphis + icols] + fElemColor[iel] * neq;
		}
	}
}

void TPZIntPointsFEM::AssembleRhsBoundary() {
	int64_t neq = fCmesh->NEquations();
	fRhsBoundary.Resize(neq, 1);
	fRhsBoundary.Zero();

	for (auto iel : fBoundaryElements) {
		TPZCompEl *cel = fCmesh->Element(iel);
		if (!cel)
			continue;
		TPZElementMatrix ef(fCmesh, TPZElementMatrix::EF);
		cel->CalcResidual(ef);
		ef.ComputeDestinationIndices();
		fRhsBoundary.AddFel(ef.fMat, ef.fSourceIndex, ef.fDestinationIndex);
	}
}

void TPZIntPointsFEM::TransferDataStructure() {

	cublasCreate(&handle_cublas);
	cusparseCreate(&handle_cusparse);

	int64_t neq = fCmesh->NEquations();
	int64_t nelem = fColSizes.size();

	cudaMalloc((void**) &dRhs, neq * sizeof(REAL));
	cudaMemset(dRhs, 0, neq * sizeof(REAL));

	cudaMalloc((void**) &dRhsBoundary, neq * sizeof(REAL));
	cudaMemcpy(dRhsBoundary, &fRhsBoundary(0, 0), neq * sizeof(REAL),
			cudaMemcpyHostToDevice);

	cudaMalloc((void**) &dSolution, neq * sizeof(REAL));
	cudaMemset(dSolution, 0, neq * sizeof(REAL));

	cudaMalloc((void**) &dPlasticStrain, fDim * fNpts * sizeof(REAL));
	cudaMemset(dPlasticStrain, 0, fDim * fNpts * sizeof(REAL));

	cudaMalloc((void**) &dStorage, fStorage.size() * sizeof(REAL));
	cudaMemcpy(dStorage, &fStorage[0], fStorage.size() * sizeof(REAL),
			cudaMemcpyHostToDevice);

	cudaMalloc((void**) &dRowSizes, nelem * sizeof(int));
	cudaMemcpy(dRowSizes, &fRowSizes[0], nelem * sizeof(int),
			cudaMemcpyHostToDevice);

	cudaMalloc((void**) &dColSizes, nelem * sizeof(int));
	cudaMemcpy(dColSizes, &fColSizes[0], nelem * sizeof(int),
			cudaMemcpyHostToDevice);

	cudaMalloc((void**) &dMatrixPosition, nelem * sizeof(int));
	cudaMemcpy(dMatrixPosition, &fMatrixPosition[0], nelem * sizeof(int),
			cudaMemcpyHostToDevice);

	cudaMalloc((void**) &dRowFirstIndex, nelem * sizeof(int));
	cudaMemcpy(dRowFirstIndex, &fRowFirstIndex[0], nelem * sizeof(int),
			cudaMemcpyHostToDevice);

	cudaMalloc((void**) &dColFirstIndex, nelem * sizeof(int));
	cudaMemcpy(dColFirstIndex, &fColFirstIndex[0], nelem * sizeof(int),
			cudaMemcpyHostToDevice);

	cudaMalloc((void**) &dIndexes, fIndexes.size() * sizeof(int));
	cudaMemcpy(dIndexes, &fIndexes[0], fIndexes.size() * sizeof(int),
			cudaMemcpyHostToDevice);

	cudaMalloc((void**) &dIndexesColor, fIndexesColor.size() * sizeof(int));
	cudaMemcpy(dIndexesColor, &fIndexesColor[0],
			fIndexesColor.size() * sizeof(int), cudaMemcpyHostToDevice);

	cudaMalloc((void**) &dWeight, fWeight.size() * sizeof(REAL));
	cudaMemcpy(dWeight, &fWeight[0], fWeight.size() * sizeof(REAL),
			cudaMemcpyHostToDevice);
}

void TPZIntPointsFEM::AssembleResidual() {
	int64_t ncolor = *std::max_element(fElemColor.begin(), fElemColor.end()) + 1;
	int64_t neq = fCmesh->NEquations();

	REAL *gather_solution;
	cudaMalloc((void**) &gather_solution, fNpts * sizeof(REAL));
	cudaMemset(gather_solution, 0, fNpts * sizeof(REAL));

	REAL *delta_strain;
	cudaMalloc((void**) &delta_strain, fDim * fNpts * sizeof(REAL));
	cudaMemset(delta_strain, 0, fDim * fNpts * sizeof(REAL));

	REAL *elastic_strain;
	cudaMalloc((void**) &elastic_strain, fDim * fNpts * sizeof(REAL));
	cudaMemset(elastic_strain, 0, fDim * fNpts * sizeof(REAL));

	REAL *sigma_trial;
	cudaMalloc((void**) &sigma_trial, fDim * fNpts * sizeof(REAL));
	cudaMemset(sigma_trial, 0, fDim * fNpts * sizeof(REAL));

	REAL *eigenvalues;
	cudaMalloc((void**) &eigenvalues, 3 * fNpts / fDim * sizeof(REAL));
	cudaMemset(eigenvalues, 0, 3 * fNpts / fDim * sizeof(REAL));

	REAL *eigenvectors;
	cudaMalloc((void**) &eigenvectors, 9 * fNpts / fDim * sizeof(REAL));
	cudaMemset(eigenvectors, 0, 9 * fNpts / fDim * sizeof(REAL));

	REAL *sigma_projected;
	cudaMalloc((void**) &sigma_projected, 3 * fNpts / fDim * sizeof(REAL));
	cudaMemset(sigma_projected, 0, 3 * fNpts / fDim * sizeof(REAL));

	REAL *sigma;
	cudaMalloc((void**) &sigma, fDim * fNpts * sizeof(REAL));
	cudaMemset(sigma, 0, fDim * fNpts * sizeof(REAL));

	REAL *nodal_forces;
	cudaMalloc((void**) &nodal_forces, fDim * fNpts * sizeof(REAL));
	cudaMemset(nodal_forces, 0, fDim * fNpts * sizeof(REAL));

	REAL *residual;
	cudaMalloc((void**) &residual, neq * ncolor * sizeof(REAL));
	cudaMemset(residual, 0, neq * ncolor * sizeof(REAL));

	cudaMemcpy(dSolution, &fSolution(0, 0), neq * sizeof(REAL), cudaMemcpyHostToDevice);
	GatherSolutionGPU(gather_solution);
	DeltaStrainGPU(gather_solution, delta_strain);
	ElasticStrainGPU(delta_strain, dPlasticStrain, elastic_strain);
	ComputeStressGPU(elastic_strain, sigma_trial);
	SpectralDecompositionGPU(sigma_trial, eigenvalues, eigenvectors);
	ProjectSigmaGPU(eigenvalues, sigma_projected);
	StressCompleteTensorGPU(sigma_projected, eigenvectors, sigma);
	NodalForcesGPU(sigma, nodal_forces);
	ColoredAssembleGPU(nodal_forces, residual);

//update strain
	ComputeStrainGPU(sigma, elastic_strain);
	PlasticStrainGPU(delta_strain, elastic_strain, dPlasticStrain);

	REAL a = 1.;
	cublasDaxpy(handle_cublas, neq, &a, &dRhsBoundary[0], 1, &residual[0], 1);

	fRhs.Resize(neq, 1);
	cudaMemcpy(&fRhs(0,0), residual, neq * sizeof(REAL), cudaMemcpyDeviceToHost);

	cudaFree(gather_solution);
	cudaFree(delta_strain);
	cudaFree(elastic_strain);
	cudaFree(sigma_trial);
	cudaFree(eigenvalues);
	cudaFree(eigenvectors);
	cudaFree(sigma_projected);
	cudaFree(sigma);
	cudaFree(nodal_forces);
	cudaFree(residual);
}

void TPZIntPointsFEM::GatherSolutionGPU(REAL *gather_solution) {
	cusparseDgthr(handle_cusparse, fDim * fNphis, dSolution, gather_solution, dIndexes, CUSPARSE_INDEX_BASE_ZERO);
}

void TPZIntPointsFEM::DeltaStrainGPU(REAL *gather_solution, REAL *delta_strain) {
	int64_t nelem = fRowSizes.size();

	int numBlocks = (nelem + NT - 1) / NT;
	DeltaStrainKernel<<<numBlocks, NT>>>(nelem, dStorage, dRowSizes, dColSizes, dMatrixPosition, dRowFirstIndex, dColFirstIndex, fNpts, fNphis, gather_solution, delta_strain);
	cudaDeviceSynchronize();
}

void TPZIntPointsFEM::ElasticStrainGPU(REAL *delta_strain, REAL *plastic_strain, REAL *elastic_strain) {
	cudaMemcpy(elastic_strain, &delta_strain[0], fDim * fNpts * sizeof(REAL), cudaMemcpyDeviceToDevice);

	REAL a = -1.;
	cublasDaxpy(handle_cublas, fDim * fNpts, &a, &plastic_strain[0], 1, &elastic_strain[0], 1);
}

void TPZIntPointsFEM::PlasticStrainGPU(REAL *delta_strain, REAL *elastic_strain, REAL *plastic_strain) {
	cudaMemcpy(plastic_strain, &delta_strain[0], fDim * fNpts * sizeof(REAL), cudaMemcpyDeviceToDevice);

	REAL a = -1.;
	cublasDaxpy(handle_cublas, fDim * fNpts, &a, &elastic_strain[0], 1, &plastic_strain[0], 1);
}

//Compute stress
void TPZIntPointsFEM::ComputeStressGPU(REAL *elastic_strain, REAL *sigma) {
	REAL lambda = fMaterial->GetPlasticModel().fER.Lambda();
	REAL mu = fMaterial->GetPlasticModel().fER.Mu();

	int numBlocks = (fNpts / fDim + NT - 1) / NT;
	ComputeStressKernel<<<numBlocks, NT>>>(fNpts, fDim, elastic_strain, sigma, mu, lambda);
	cudaDeviceSynchronize();
}

//Compute strain
void TPZIntPointsFEM::ComputeStrainGPU(REAL *sigma, REAL *elastic_strain) {
	REAL E = fMaterial->GetPlasticModel().fER.E();
	REAL nu = fMaterial->GetPlasticModel().fER.Poisson();

	int numBlocks = (fNpts / fDim + NT - 1) / NT;
	ComputeStrainKernel<<<numBlocks, NT>>>(fNpts, fDim, sigma, elastic_strain, nu, E, dWeight);
	cudaDeviceSynchronize();
}

void TPZIntPointsFEM::SpectralDecompositionGPU(REAL *sigma_trial, REAL *eigenvalues, REAL *eigenvectors) {
	int numBlocks = (fNpts / fDim + NT - 1) / NT;
	SpectralDecompositionKernel<<<numBlocks, NT>>>(fNpts, fDim, sigma_trial, eigenvalues, eigenvectors);
	cudaDeviceSynchronize();
}

void TPZIntPointsFEM::ProjectSigmaGPU(REAL *eigenvalues, REAL *sigma_projected) {

	REAL mc_phi = fMaterial->GetPlasticModel().fYC.Phi();
	REAL mc_psi = fMaterial->GetPlasticModel().fYC.Psi();
	REAL mc_cohesion = fMaterial->GetPlasticModel().fYC.Cohesion();
	REAL K = fMaterial->GetPlasticModel().fER.K();
	REAL G = fMaterial->GetPlasticModel().fER.G();

	REAL *m_type;
	cudaMalloc((void**) &m_type, fNpts / fDim * sizeof(REAL));
	cudaMemset(m_type, 0, fNpts / fDim * sizeof(REAL));

	REAL *alpha;
	cudaMalloc((void**) &alpha, fNpts / fDim * sizeof(REAL));
	cudaMemset(alpha, 0, fNpts / fDim * sizeof(REAL));

	int numBlocks = (fNpts / fDim + NT - 1) / NT;
	ProjectSigmaKernel<<<numBlocks, NT>>>(fNpts, fDim, mc_phi, mc_psi, mc_cohesion, K, G, eigenvalues, sigma_projected, m_type, alpha);
	cudaDeviceSynchronize();

}

void TPZIntPointsFEM::StressCompleteTensorGPU(REAL *sigma_projected, REAL *eigenvectors, REAL *sigma) {
	int numBlocks = (fNpts / fDim + NT - 1) / NT;
	StressCompleteTensorKernel<<<numBlocks, NT>>>(fNpts, fDim, sigma_projected, eigenvectors, sigma, dWeight);
	cudaDeviceSynchronize();


}

void TPZIntPointsFEM::NodalForcesGPU(REAL *sigma, REAL *nodal_forces) {
	int64_t nelem = fRowSizes.size();

	int numBlocks = (nelem + NT - 1) / NT;
	NodalForcesKernel<<<numBlocks, NT>>>(nelem, dStorage, dRowSizes, dColSizes, dMatrixPosition, dRowFirstIndex, dColFirstIndex, fNpts, fNphis, sigma, nodal_forces);
	cudaDeviceSynchronize();
}

void TPZIntPointsFEM::ColoredAssembleGPU(REAL *nodal_forces, REAL *residual) {
	int64_t ncolor = *std::max_element(fElemColor.begin(), fElemColor.end())
			+ 1;
	int64_t sz = fIndexes.size();
	int64_t neq = fCmesh->NEquations();

	cusparseDsctr(handle_cusparse, sz, &nodal_forces[0], &dIndexesColor[0], &residual[0], CUSPARSE_INDEX_BASE_ZERO);

	int64_t colorassemb = ncolor / 2.;
	REAL alpha = 1.;
	while (colorassemb > 0) {

		int64_t firsteq = (ncolor - colorassemb) * neq;
		cublasDaxpy(handle_cublas, firsteq, &alpha, &residual[firsteq], 1., &residual[0], 1.);

		ncolor -= colorassemb;
		colorassemb = ncolor / 2;
	}
}
