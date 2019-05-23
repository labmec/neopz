#include "pzreal.h"

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
