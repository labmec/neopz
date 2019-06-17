#include "pzreal.h"

__device__ void NormalizeDevice(REAL *sigma, REAL &maxel) {
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

__device__ void IntervalDevice(REAL *sigma, REAL *interval) {
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

__device__ void NewtonIterationsDevice(REAL *interval, REAL *sigma, REAL *eigenvalues, REAL &maxel) {
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

__device__ void Multiplicity1Device(REAL *sigma, REAL eigenvalue, REAL *eigenvector) {
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

__device__ void Multiplicity2Device(REAL *sigma, REAL eigenvalue, REAL *eigenvector1,
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

__device__ void EigenvectorsDevice(REAL *sigma, REAL *eigenvalues, REAL *eigenvectors,
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
			Multiplicity1Device(sigma, eigenvalues[0], &eigenvectors[0]);
		} else if (eigenvalues[0] == eigenvalues[1]) {
			Multiplicity2Device(sigma, eigenvalues[0], &eigenvectors[0], &eigenvectors[3]);
		} else if (eigenvalues[0] == eigenvalues[2]) {
			Multiplicity2Device(sigma, eigenvalues[0], &eigenvectors[0], &eigenvectors[6]);
		}
		if (eigenvalues[1] != eigenvalues[0] && eigenvalues[1] != eigenvalues[2]) {
			Multiplicity1Device(sigma, eigenvalues[1], &eigenvectors[3]);
		} else if (eigenvalues[1] == eigenvalues[2]) {
			Multiplicity2Device(sigma, eigenvalues[1], &eigenvectors[3], &eigenvectors[6]);
		}
		if (eigenvalues[2] != eigenvalues[0] && eigenvalues[2] != eigenvalues[1]) {
			Multiplicity1Device(sigma, eigenvalues[2], &eigenvectors[6]);
		}
	}
}

__global__ void SpectralDecompositionKernel(REAL *sigma_trial, REAL *eigenvalues, REAL *eigenvectors, int64_t npts) {
	int ipts = blockIdx.x * blockDim.x + threadIdx.x;

	__shared__ REAL maxel;
	__shared__ REAL interval[2];
	if (ipts < npts) {
		NormalizeDevice(&sigma_trial[4 * ipts], maxel);
		IntervalDevice(&sigma_trial[4 * ipts], &interval[0]);
		NewtonIterationsDevice(&interval[0], &sigma_trial[4 * ipts], &eigenvalues[3 * ipts], maxel);
		EigenvectorsDevice(&sigma_trial[4 * ipts], &eigenvalues[3 * ipts], &eigenvectors[9 * ipts], maxel);
	}
}
