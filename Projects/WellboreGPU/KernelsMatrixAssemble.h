#include "pzreal.h"

#ifdef O_LINEAR
#define ndof 8
#define NT_sm 256
#elif O_QUADRATIC
#define ndof 18
#define NT_sm 110
#elif O_CUBIC
#define ndof 32
#define NT_sm 64
#endif

__constant__ REAL De[3 * 3];

__device__ void ComputeTangentMatrixDevice(int el_npts, int el_dofs, REAL *storage, REAL weight, REAL *K){   

    REAL DeBip[3 * ndof];
    for(int i = 0; i < 3 * el_dofs; i++) DeBip[i] = 0.0;
    REAL omega = weight;
    MultAddDevice(false, 3, el_dofs, 3, De, storage, DeBip, 1., 0.);
    MultAddDevice(true, el_dofs, el_dofs, 3, storage, DeBip, K, omega, 1.);
}


__device__ int64_t me(int *ia_to_sequence, int *ja_to_sequence, int64_t & i_dest, int64_t & j_dest) {
    int64_t row(i_dest),col(j_dest);
    if (i_dest > j_dest) {
        int64_t temp = i_dest;
        row = col;
        col = temp;
    }
    for(int ic=ia_to_sequence[row] ; ic < ia_to_sequence[row+1]; ic++ ) {
        if ( ja_to_sequence[ic] == col )
        {
            return ic;
        }
    }
    return 0; 
}

__global__ 
void MatrixAssembleKernel(int nel, REAL *Kg, int64_t *el_color_index, REAL *weight, int *dof_indexes, 
    REAL *storage, int *rowsizes, int *colsizes, int *rowfirstindex, int *colfirstindex, int *matrixposition, int *ia_to_sequence, int *ja_to_sequence) {

    int tid = blockIdx.x * blockDim.x + threadIdx.x;

    if(tid < nel) {
        int iel = el_color_index[tid];

        int el_npts = rowsizes[iel]/3;
        int el_dofs = colsizes[iel];
        int colpos = colfirstindex[iel];
        int first_el_ip = rowfirstindex[iel]/3;
        int matpos = matrixposition[iel];

        int64_t dest[ndof];
        REAL K[ndof * ndof];
        for(int i = 0; i < el_dofs; i++) dest[i] = dof_indexes[colpos + i];
        for(int i = 0; i < el_dofs * el_dofs; i++) K[i] = 0;
        __shared__ REAL s_storage[NT_sm * ndof * 3]; // max allowed word is 48k
        for (int ip = 0; ip < el_npts; ip++) {
            for(int i = 0; i < ndof * 3; i++) {
               s_storage[i + threadIdx.x * ndof * 3] = storage[matpos + i + ip * ndof * 3];
            }
            ComputeTangentMatrixDevice(el_npts, el_dofs, &s_storage[threadIdx.x * ndof * 3], weight[first_el_ip + ip], K);   
        }
 
        for (int i_dof = 0; i_dof < el_dofs; i_dof++) {
            int64_t i_dest = dest[i_dof];
            for (int j_dof = i_dof; j_dof < el_dofs; j_dof++) {
                int64_t j_dest = dest[j_dof];
                int64_t index = me(ia_to_sequence, ja_to_sequence, i_dest, j_dest);
                Kg[index] += K[i_dof * el_dofs + j_dof];
            }
        }           
     }
}


__global__ 
void MatrixAssembleKernelGS(int nel, REAL *Kc, int64_t *el_color_index, REAL *weight, int *dof_indexes, 
	REAL *storage, int *rowsizes, int *colsizes, int *rowfirstindex, int *colfirstindex, int *matrixposition) {

	// int tid = blockIdx.x;
    int tid = blockIdx.x * blockDim.x + threadIdx.x;

    if(tid < nel) {
        int iel = el_color_index[tid];

        int el_npts = rowsizes[iel]/3;
        int el_dofs = colsizes[iel];
        int colpos = colfirstindex[iel];
        int first_el_ip = rowfirstindex[iel]/3;
        int matpos = matrixposition[iel];

        REAL K[ndof * ndof];
        for(int i = 0; i < el_dofs * el_dofs; i++) K[i] = 0;
        __shared__ REAL s_storage[NT_sm * ndof * 3]; // max allowed word is 48k
        for (int ip = 0; ip < el_npts; ip++) {
            for(int i = 0; i < ndof * 3; i++) {
               s_storage[i + threadIdx.x * ndof * 3] = storage[matpos + i + ip * ndof * 3];
            }
            ComputeTangentMatrixDevice(el_npts, el_dofs, &s_storage[threadIdx.x * ndof * 3], weight[first_el_ip + ip], K);   
        }

        int stride = tid*(el_dofs * el_dofs + el_dofs)/2;
        int c = stride;
        for(int i_dof = 0 ; i_dof < el_dofs; i_dof++){
            for(int j_dof = i_dof; j_dof < el_dofs; j_dof++){
                Kc[c] += K[i_dof * el_dofs + j_dof];
                c++;
            }
        }	
     }
}

