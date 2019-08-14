#include "pzreal.h"

#ifdef __CUDACC__
__device__ 
#endif
bool PhiPlane(REAL *eigenvalues, REAL *sigma_projected, REAL mc_phi, REAL mc_cohesion) {
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

#ifdef __CUDACC__
__device__ 
#endif
bool ReturnMappingMainPlane(REAL *sigma_trial, REAL *sigma_projected, REAL &m_hardening, REAL mc_phi, REAL mc_psi, REAL mc_cohesion, REAL K, REAL G) {
    REAL eigenvalues[3];
    eigenvalues[0] = sigma_trial[0];
    eigenvalues[1] = sigma_trial[1];
    eigenvalues[2] = sigma_trial[2];

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

#ifdef __CUDACC__
__device__
#endif
void ComputePlaneTangent(REAL *tang, REAL mc_phi, REAL mc_psi, REAL K, REAL G) {
//    for (int i = 0; i < 9; i++) tang[i] = 0.;
//    tang[0] = 1.;
//    tang[4] = 1.;
//    tang[8] = 1.;
//    return;
	const REAL sin_phi = sin(mc_phi);
	const REAL sin_psi = sin(mc_psi);
	const REAL denominator = 6.0 * G + 2.0 * (G + 3.0 * K) * sin_phi * sin_psi;

	// First column
	tang[0] = (sin_phi - 1.0) * (-3.0 * G + (G + 3.0 * K) * sin_psi) / denominator;
	tang[3] = (2.0 * G - 3.0 * K) * (sin_phi + 1.0) * sin_psi / denominator;
	tang[6] = -(sin_phi + 1.0) * (-3.0 * G + (G + 3.0 * K) * sin_psi) / denominator;

	// Second column
	tang[1] = 0.0;
	tang[4] = 1.0;
	tang[7] = 0.0;

	// Third column
	tang[2] = -(sin_phi - 1.0) * (3.0 * G + (G + 3.0 * K) * sin_psi) / denominator;
	tang[5] = (2.0 * G - 3.0 * K) * (sin_phi - 1.0) * sin_psi / denominator;
	tang[8] = (sin_phi + 1.0) * (3.0 * G + (G + 3.0 * K) * sin_psi) / denominator;

}

#ifdef __CUDACC__
__device__ 
#endif
bool ReturnMappingRightEdge(REAL *sigma_trial, REAL *sigma_projected, REAL &m_hardening, REAL mc_phi, REAL mc_psi, REAL mc_cohesion, REAL K, REAL G) {
    REAL eigenvalues[3];
    eigenvalues[0] = sigma_trial[0];
    eigenvalues[1] = sigma_trial[1];
    eigenvalues[2] = sigma_trial[2];

    const REAL sinphi = sin(mc_phi);
	const REAL sinpsi = sin(mc_psi);
	const REAL cosphi = cos(mc_phi);

	 REAL gamma[2], phi[2], sigma_bar[2], ab[2];

	 REAL jac[2][2], jac_inv[2][2];

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

#ifdef __CUDACC__
__device__
#endif
void ComputeRightEdgeTangent(REAL *tang, REAL mc_phi, REAL mc_psi, REAL K, REAL G) {
//    for (int i = 0; i < 9; i++) tang[i] = 0.;
//    tang[0] = 1.;
//    tang[4] = 1.;
//    tang[8] = 1.;
//    return;
	const REAL sin_phi = sin(mc_phi);
	const REAL sin_psi = sin(mc_psi);
	const REAL a = 4.0 * G * (1.0 + (1.0/3.0) * sin_phi * sin_psi) + 4.0 * K * sin_phi * sin_psi;
	const REAL b = 2.0 * G * (1.0 + sin_phi + sin_psi - (1.0/3.0) * sin_phi * sin_psi) + 4.0 * K * sin_phi * sin_psi;

	// First column
	tang[0] = (3.0*a + 3.0*b - 4.0*(1.0 + sin_phi)*(3.0*K*sin_psi + G*(3.0 + sin_psi)))/(3.*(a + b));
	tang[3] = (2.0*(1.0 + sin_phi)*(-6.0*K*sin_psi + G*(3.0 + sin_psi)))/(3.*(a + b));
	tang[6] = (2.0*(1.0 + sin_phi)*(-6.0*K*sin_psi + G*(3.0 + sin_psi)))/(3.*(a + b));

	// Second column
	tang[1] = (-2*(-1 + sin_phi)*(3*K*sin_psi + G*(3 + sin_psi)))/(3.*(a + b));
	tang[4] = (-3*b*b + 3*a*(a + 2*G*(-1 + sin_phi)) -
				  2*(a*G + 2*b*G + 3*a*K - 3*b*K)*(-1 + sin_phi)*sin_psi)/(3.*(a - b)*(a + b));
	tang[7] = (2*(-1 + sin_phi)*(b*G*(-3 + sin_psi) + a*(2*G - 3*K)*sin_psi + 3*b*K*sin_psi))/
				 (3.*(a*a - b*b));

	// Third column
	tang[2] = (-2*(-1 + sin_phi)*(3*K*sin_psi + G*(3 + sin_psi)))/(3.*(a + b));
	tang[5] = (2*(-1 + sin_phi)*(b*G*(-3 + sin_psi) + a*(2*G - 3*K)*sin_psi + 3*b*K*sin_psi))/
				 (3.*(a*a - b*b));
	tang[8] = (-3*b*b + 3*a*(a + 2*G*(-1 + sin_phi)) -
				  2*(a*G + 2*b*G + 3*a*K - 3*b*K)*(-1 + sin_phi)*sin_psi)/(3.*(a - b)*(a + b));

}

#ifdef __CUDACC__
__device__ 
#endif
bool ReturnMappingLeftEdge(REAL *sigma_trial, REAL *sigma_projected, REAL &m_hardening, REAL mc_phi, REAL mc_psi, REAL mc_cohesion, REAL K, REAL G) {
    REAL eigenvalues[3];
    eigenvalues[0] = sigma_trial[0];
    eigenvalues[1] = sigma_trial[1];
    eigenvalues[2] = sigma_trial[2];

    const REAL sinphi = sin(mc_phi);
	const REAL sinpsi = sin(mc_psi);
	const REAL cosphi = cos(mc_phi);
	const REAL sinphi2 = sinphi * sinphi;
	const REAL cosphi2 = 1. - sinphi2;

	 REAL gamma[2], phi[2], sigma_bar[2], ab[2];

	 REAL jac[2][2], jac_inv[2][2];

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

#ifdef __CUDACC__
__device__
#endif
void ComputeLeftEdgeTangent(REAL *tang, REAL mc_phi, REAL mc_psi, REAL K, REAL G) {
//    for (int i = 0; i < 9; i++) tang[i] = 0.;
//    tang[0] = 1.;
//    tang[4] = 1.;
//    tang[8] = 1.;
//    return;
	const REAL sin_phi = sin(mc_phi);
	const REAL sin_psi = sin(mc_psi);
	const REAL a = 4.0 * G * (1.0 + (1.0/3.0) * sin_phi * sin_psi) + 4.0 * K * sin_phi * sin_psi;
	const REAL b = 2.0 * G * (1.0 - sin_phi - sin_psi - (1.0/3.0) * sin_phi * sin_psi) + 4.0 * K * sin_phi * sin_psi;

	// First column
	tang[0] = (-3*b*b + 3*a*(a - 2*G*(1 + sin_phi)) -
				  2*(a*G + 2*b*G + 3*a*K - 3*b*K)*(1 + sin_phi)*sin_psi)/(3.*(a - b)*(a + b));
	tang[3] = (2*(1 + sin_phi)*(a*(2*G - 3*K)*sin_psi + 3*b*K*sin_psi + b*G*(3 + sin_psi)))/
				 (3.*(a*a - b*b));
	tang[6] = (-2*(1 + sin_phi)*(G*(-3 + sin_psi) + 3*K*sin_psi))/(3.*(a + b));

	// Second column
	tang[1] = (2*(1 + sin_phi)*(a*(2*G - 3*K)*sin_psi + 3*b*K*sin_psi + b*G*(3 + sin_psi)))/
				 (3.*(a*a - b*b));
	tang[4] = (-3*b*b + 3*a*(a - 2*G*(1 + sin_phi)) -
				  2*(a*G + 2*b*G + 3*a*K - 3*b*K)*(1 + sin_phi)*sin_psi)/(3.*(a - b)*(a + b));
	tang[7] = (-2*(1 + sin_phi)*(G*(-3 + sin_psi) + 3*K*sin_psi))/(3.*(a + b));

	// Third column
	tang[2] = (2*(-1 + sin_phi)*(G*(-3 + sin_psi) - 6*K*sin_psi))/(3.*(a + b));
	tang[5] = (2*(-1 + sin_phi)*(G*(-3 + sin_psi) - 6*K*sin_psi))/(3.*(a + b));
	tang[8] = (3*a + 3*b - 4*(-1 + sin_phi)*(G*(-3 + sin_psi) + 3*K*sin_psi))/(3.*(a + b));

}

#ifdef __CUDACC__
__device__ 
#endif
void ReturnMappingApex(REAL *sigma_trial, REAL *sigma_projected, REAL &m_hardening, REAL mc_phi, REAL mc_psi, REAL mc_cohesion, REAL K) {
	const REAL cotphi = 1. / tan(mc_phi);

	REAL ptrnp1 = 0.;
	for (int i = 0; i < 3; i++) {
		ptrnp1 += sigma_trial[i];
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

#ifdef __CUDACC__
__device__
#endif
void ComputeApexGradient(REAL *gradient) {
//    for (int i = 0; i < 9; i++) gradient[i] = 0.;
//    gradient[0] = 1.;
//    gradient[4] = 1.;
//    gradient[8] = 1.;
//    return;
	for (int i = 0; i < 3 * 3; i++) {
		gradient[i] = 0;
	}
}

#ifdef __CUDACC__
__device__ 
#endif
void ProjectSigma(REAL *eigenvalues, REAL *sigma_projected, int &m_type, REAL &alpha, REAL mc_phi, REAL mc_psi, REAL mc_cohesion, REAL K, REAL G, REAL *gradient, bool require_gradient_Q) {

	bool check_Q = false;

	m_type = 0;
	check_Q = PhiPlane(eigenvalues, sigma_projected, mc_phi, mc_cohesion); //elastic domain
	if(check_Q) {
		if (require_gradient_Q) {
			for (int i = 0; i < 9; i++) gradient[i] = 0.;
			gradient[0] = 1.;
			gradient[4] = 1.;
			gradient[8] = 1.;
		}
		return;
	}
	else { //plastic domain
		m_type = 1;
		check_Q = ReturnMappingMainPlane(eigenvalues, sigma_projected, alpha, mc_phi, mc_psi, mc_cohesion, K, G); //main plane
		if (check_Q && require_gradient_Q) {
			ComputePlaneTangent(gradient, mc_phi, mc_psi, K, G);
		}
		if (!check_Q) { //edges or apex
			if (((1 - sin(mc_psi)) * eigenvalues[0] - 2. * eigenvalues[1] + (1 + sin(mc_psi)) * eigenvalues[2]) > 0) { // right edge
				check_Q = ReturnMappingRightEdge(eigenvalues, sigma_projected, alpha, mc_phi, mc_psi, mc_cohesion, K, G);
				if (check_Q && require_gradient_Q) {
					ComputeRightEdgeTangent(gradient, mc_phi, mc_psi, K, G);
				}
			} else { //left edge
				check_Q = ReturnMappingLeftEdge(eigenvalues, sigma_projected, alpha, mc_phi, mc_psi, mc_cohesion, K, G);
				if (check_Q && require_gradient_Q) {
					ComputeLeftEdgeTangent(gradient, mc_phi, mc_psi, K, G);
				}
			}
			if (!check_Q) { //apex
				m_type = -1;
				ReturnMappingApex(eigenvalues, sigma_projected, alpha, mc_phi, mc_psi, mc_cohesion, K);
				if (require_gradient_Q) {
					ComputeApexGradient(gradient);
				}
			}
		}
	}
}
