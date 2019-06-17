#include "pzreal.h"

__global__ void ComputeStressKernel(REAL *elastic_strain, REAL *sigma, int64_t npts, REAL mu, REAL lambda) {
	int ipts = blockIdx.x * blockDim.x + threadIdx.x;

//	extern __shared__ REAL elastic_strain_s[];
//	memcpy(&elastic_strain_s[0], &elastic_strain[2 * blockIdx.x * blockDim.x], 2 * blockDim.x * sizeof(REAL));
//	memcpy(&elastic_strain_s[2 * blockDim.x], &elastic_strain[fNpts + 2 * blockIdx.x * blockDim.x], 2 * blockDim.x * sizeof(REAL));
//	__syncthreads();
//	if (ipts < blockDim.x) {
//		sigma[4 * ipts] = elastic_strain_s[2 * ipts] * (lambda + 2. * mu) + elastic_strain[2 * ipts + 2 * blockDim.x + 1] * lambda; // Sigma xx
//		sigma[4 * ipts + 1] = elastic_strain[2 * ipts + 2 * blockDim.x + 1] * (lambda + 2. * mu) + elastic_strain_s[2 * ipts] * lambda; // Sigma yy
//		sigma[4 * ipts + 2] = lambda * (elastic_strain_s[2 * ipts] + elastic_strain[2 * ipts + 2 * blockDim.x + 1]); // Sigma zz
//		sigma[4 * ipts + 3] = mu * (elastic_strain_s[2 * ipts + 1] + elastic_strain[2 * ipts + 2 * blockDim.x]); // Sigma xy
//	}

	if (ipts < npts) {
		sigma[4 * ipts] = elastic_strain[2 * ipts] * (lambda + 2. * mu) + elastic_strain[2 * ipts + 2 * npts + 1] * lambda; // Sigma xx
		sigma[4 * ipts + 1] = elastic_strain[2 * ipts + 2 * npts + 1] * (lambda + 2. * mu) + elastic_strain[2 * ipts] * lambda; // Sigma yy
		sigma[4 * ipts + 2] = lambda * (elastic_strain[2 * ipts] + elastic_strain[2 * ipts + 2 * npts + 1]); // Sigma zz
		sigma[4 * ipts + 3] = mu * (elastic_strain[2 * ipts + 1] + elastic_strain[2 * ipts + 2 * npts]); // Sigma xy
	}
}

// __global__ void ComputeStrainKernel(int64_t fNpts, int fDim, REAL *sigma, REAL *elastic_strain, REAL nu, REAL E, REAL *weight) {
// 	int ipts = blockIdx.x * blockDim.x + threadIdx.x;

// 	if (ipts < fNpts / fDim) {
// 		elastic_strain[2 * ipts + 0] = 1 / weight[ipts] * (1. / E * (sigma[2 * ipts] * (1. - nu * nu) - sigma[2 * ipts + fNpts + 1] * (nu + nu * nu))); //exx
// 		elastic_strain[2 * ipts + 1] = 1 / weight[ipts] * ((1. + nu) / E * sigma[2 * ipts + 1]); //exy
// 		elastic_strain[2 * ipts + fNpts + 0] = elastic_strain[2 * ipts + 1]; //exy
// 		elastic_strain[2 * ipts + fNpts + 1] = 1 / weight[ipts] * (1. / E * (sigma[2 * ipts + fNpts + 1] * (1. - nu * nu) - sigma[2 * ipts] * (nu + nu * nu))); //eyy
// 	}
// }

// __global__ void StressCompleteTensorKernel(int64_t fNpts, int fDim, REAL *sigma_projected, REAL *eigenvectors, REAL *sigma, REAL *weight) {
// 	int ipts = blockIdx.x * blockDim.x + threadIdx.x;

// 	if (ipts < fNpts / fDim) {

// 		sigma[2 * ipts + 0] = weight[ipts] * (sigma_projected[3 * ipts + 0] * eigenvectors[9 * ipts + 0] * eigenvectors[9 * ipts + 0] + sigma_projected[3 * ipts + 1]
// 								* eigenvectors[9 * ipts + 3] * eigenvectors[9 * ipts + 3] + sigma_projected[3 * ipts + 2] * eigenvectors[9 * ipts + 6] * eigenvectors[9 * ipts + 6]);
// 		sigma[2 * ipts + 1] = weight[ipts] * (sigma_projected[3 * ipts + 0] * eigenvectors[9 * ipts + 0] * eigenvectors[9 * ipts + 1] + sigma_projected[3 * ipts + 1]
// 								* eigenvectors[9 * ipts + 3] * eigenvectors[9 * ipts + 4] + sigma_projected[3 * ipts + 2] * eigenvectors[9 * ipts + 6] * eigenvectors[9 * ipts + 7]);
// 		sigma[2 * ipts + fNpts] = sigma[2 * ipts + 1];
// 		sigma[2 * ipts + fNpts + 1] = weight[ipts] * (sigma_projected[3 * ipts + 0] * eigenvectors[9 * ipts + 1] * eigenvectors[9 * ipts + 1] + sigma_projected[3 * ipts + 1]
// 								* eigenvectors[9 * ipts + 4] * eigenvectors[9 * ipts + 4] + sigma_projected[3 * ipts + 2] * eigenvectors[9 * ipts + 7] * eigenvectors[9 * ipts + 7]);
// 	}
// }

