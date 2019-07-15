#include "pzreal.h"
#include "TPZTensor.h"

#ifdef __CUDACC__
__device__ 
#endif
void Normalize(REAL *sigma, REAL &maxel) {
	maxel = sigma[0];
	for (int i = 1; i < 6; i++) {
		if (fabs(sigma[i]) > fabs(maxel)) {
			maxel = sigma[i];
		}
	}
	if(maxel != 0) {
		for (int i = 0; i < 6; i++) {
			sigma[i] /= maxel;
		}
	}

}

#ifdef __CUDACC__
__device__ 
#endif
void Interval(REAL *sigma, REAL *interval) {
	 REAL lower_vec[3];
	 REAL upper_vec[3];

	//row 1 |sigma_xx sigma_xy 0|
	REAL row_1 = sigma[_XY_] + sigma[_XZ_];
    lower_vec[0] = sigma[_XX_] - fabs(row_1);
    upper_vec[0] = sigma[_XX_] + fabs(row_1);

	//row 2 |sigma_xy sigma_yy 0|
    REAL row_2 = sigma[_XY_] + sigma[_YZ_];
    lower_vec[1] = sigma[_YY_] - fabs(row_2);
    upper_vec[1] = sigma[_YY_] + fabs(row_2);

	//row 3 |0 0 sigma_zz|
    REAL row_3 = sigma[_XY_] + sigma[_YZ_];
    lower_vec[2] = sigma[_ZZ_] - fabs(row_3);
    upper_vec[2] = sigma[_ZZ_] + fabs(row_3);

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

#ifdef __CUDACC__
__device__ 
#endif
void NewtonIterations(REAL *interval, REAL *sigma, REAL *eigenvalues, REAL &maxel) {
	int numiterations = 20;
	REAL tol = 10e-12;

	REAL res, f, df, x;
	int it;

	for (int i = 0; i < 2; i++) {
		x = interval[i];
		it = 0;

		REAL sigmaxy2 = sigma[_XY_]*sigma[_XY_];
        REAL sigmaxz2 = sigma[_XZ_]*sigma[_XZ_];
        REAL sigmayz2 = sigma[_YZ_]*sigma[_YZ_];

		f = -x*x*x - sigmaxz2* sigma[_YY_] - sigmaxy2*sigma[_ZZ_] - sigmayz2*sigma[_XX_] + 2* sigma[_XY_]*sigma[_XZ_]*sigma[_YZ_] +
        	sigma[_XX_]*sigma[_YY_]*sigma[_ZZ_] + x*x*(sigma[_XX_] + sigma[_YY_] + sigma[_ZZ_]) + x*(sigmaxy2 + sigmaxz2 + sigmayz2 -
        	sigma[_XX_]*sigma[_YY_] - sigma[_XX_]*sigma[_ZZ_] - sigma[_YY_]*sigma[_ZZ_]);

		res = abs(f);

		while (it < numiterations && res > tol) {
            df = -3*x*x + sigmaxy2 + sigmaxz2 + sigmayz2 - sigma[_XX_]*sigma[_YY_] - sigma[_XX_]*sigma[_ZZ_] - sigma[_YY_]*sigma[_ZZ_] + 2*x*(sigma[_XX_] + sigma[_YY_] + sigma[_ZZ_]);

			x -= f / df;
            f = -x*x*x - sigmaxz2* sigma[_YY_] - sigmaxy2*sigma[_ZZ_] - sigmayz2*sigma[_XX_] + 2* sigma[_XY_]*sigma[_XZ_]*sigma[_YZ_] +
                sigma[_XX_]*sigma[_YY_]*sigma[_ZZ_] + x*x*(sigma[_XX_] + sigma[_YY_] + sigma[_ZZ_]) + x*(sigmaxy2 + sigmaxz2 + sigmayz2 -
                sigma[_XX_]*sigma[_YY_] - sigma[_XX_]*sigma[_ZZ_] - sigma[_YY_]*sigma[_ZZ_]);
			res = abs(f);
			it++;
		}
		eigenvalues[i] = x;

	}
	eigenvalues[2] = sigma[_XX_] + sigma[_YY_] + sigma[_ZZ_] - eigenvalues[0] - eigenvalues[1];

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

#ifdef __CUDACC__
__device__ 
#endif
void Multiplicity1(REAL *sigma, REAL eigenvalue, REAL *eigenvector) {
	 REAL det[3];
    det[0] = (sigma[_XX_] - eigenvalue) * (sigma[_YY_] - eigenvalue) - sigma[_XY_] * sigma[_XY_];
    det[1] = (sigma[_XX_] - eigenvalue) * (sigma[_ZZ_] - eigenvalue) - sigma[_XZ_] * sigma[_XZ_];
    det[2] = (sigma[_YY_] - eigenvalue) * (sigma[_ZZ_] - eigenvalue) - sigma[_YZ_] * sigma[_YZ_];

	REAL maxdet = fabs(det[0]);
	for (int i = 1; i < 3; i++) {
		if (fabs(det[i]) > fabs(maxdet)) {
			maxdet = fabs(det[i]);
		}
	}
	 REAL v[3];
    if (maxdet == fabs(det[0])) {
        v[0] = 1 / det[0] * (-(sigma[_YY_] - eigenvalue) * sigma[_XZ_] + sigma[_XY_] * sigma[_YZ_]);
        v[1] = 1 / det[0] * (-(sigma[_XX_] - eigenvalue) * sigma[_YZ_] + sigma[_XY_] * sigma[_XZ_]);
        v[2] = 1;

    }
    else if (maxdet == fabs(det[1])) {
        v[0] = 1 / det[1] * (-(sigma[_ZZ_] - eigenvalue) * sigma[_XY_] + sigma[_XZ_] * sigma[_YZ_]);
        v[1] = 1;
        v[2] = 1 / det[1] * (-(sigma[_XX_] - eigenvalue) * sigma[_YZ_] + sigma[_XY_] * sigma[_XZ_]);

    }
    else {
        v[0] = 1;
        v[1] = 1 / det[2] * (-(sigma[_ZZ_] - eigenvalue) * sigma[_XY_] + sigma[_XZ_] * sigma[_YZ_]);
        v[2] = 1 / det[2] * (-(sigma[_YY_] - eigenvalue) * sigma[_XZ_] + sigma[_XY_] * sigma[_YZ_]);
    }
	REAL norm = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
	eigenvector[0] = v[0] / norm;
	eigenvector[1] = v[1] / norm;
	eigenvector[2] = v[2] / norm;
}

#ifdef __CUDACC__
__device__ 
#endif
void Multiplicity2(REAL *sigma, REAL eigenvalue, REAL *eigenvector1,
		REAL *eigenvector2) {
	 REAL x[3];
    x[0] = sigma[_XX_] - eigenvalue;
    x[1] = sigma[_YY_] - eigenvalue;
    x[2] = sigma[_ZZ_] - eigenvalue;

	REAL maxx = fabs(x[0]);
	for (int i = 1; i < 3; i++) {
		if (fabs(x[i]) > fabs(maxx)) {
			maxx = fabs(x[i]);
		}
	}

	 REAL v1[3];
	 REAL v2[3];

    if (maxx == fabs(x[0])) {
        v1[0] = -sigma[_XY_] / x[0];
        v1[1] = 1;
        v1[2] = 0;

        v2[0] = -sigma[_XZ_] / x[0];
        v2[1] = 0;
        v2[2] = 1;

    }
    else if (maxx == fabs(x[1])) {
        v1[0] = 1;
        v1[1] = -sigma[_XY_] / x[1];
        v1[2] = 0;

        v2[0] = 0;
        v2[1] = -sigma[_YZ_] / x[1];
        v2[2] = 1;

    }
    else {
        v1[0] = 1;
        v1[1] = 0;
        v1[2] = -sigma[_XZ_] / x[2];

        v2[0] = 0;
        v2[1] = 1;
        v2[2] = -sigma[_YZ_] / x[2];

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

#ifdef __CUDACC__
__device__ 
#endif
void Eigenvectors(REAL *sigma, REAL *eigenvalues, REAL *eigenvectors,
		REAL &maxel) {
    sigma[_XX_]*=maxel;
    sigma[_YY_]*=maxel;
    sigma[_ZZ_]*=maxel;
    sigma[_XY_]*=maxel;
    sigma[_XZ_]*=maxel;
    sigma[_YZ_]*=maxel;

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

#ifdef __CUDACC__
__device__ 
#endif
void SpectralDecomposition(REAL *sigma_trial, REAL *eigenvalues, REAL *eigenvectors) {
	REAL maxel;
	REAL interval[2];
	Normalize(sigma_trial, maxel);
	Interval(sigma_trial, interval);
	NewtonIterations(interval, sigma_trial, eigenvalues, maxel);
	Eigenvectors(sigma_trial, eigenvalues, eigenvectors, maxel);
}
