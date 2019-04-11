/**
 * @file
 * @brief Contains the TPZSolveMatrix class which implements a solution based on a matrix procedure.
 */

#ifndef TPZSolveMatrix_h
#define TPZSolveMatrix_h
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzinterpolationspace.h"
#include "pzcmesh.h"
#include "TElastoPlasticData.h"

#ifdef USING_MKL
#include "mkl.h"
#endif

#ifdef __CUDACC__
#include <cuda.h>
#include <cublas_v2.h>
#include <cusparse.h>
#endif

class TPZSolveMatrix {

public:

    TPZSolveMatrix() {
        fTotalStrain.Resize(0,0);
        fPlasticStrain.Resize(0,0);
        fSolution.Resize(0,0);
        fNpts = -1;
        fNphis = -1;
        fStorage.resize(0);
        fColSizes.resize(0);
        fRowSizes.resize(0);
        fMatrixPosition.resize(0);
        fColFirstIndex.resize(0);
        fRowFirstIndex.resize(0);
        fElemColor.resize(0);
        fIndexes.resize(0);
        fIndexesColor.resize(0);
        fWeight.resize(0);

    }

    TPZSolveMatrix(TPZCompMesh *cmesh, TElastoPlasticData materialdata) {
        SetCompMesh(cmesh);
        SetMaterialData(materialdata);
        SetDataStructure();
    }

    ~TPZSolveMatrix() {

    }

    TPZSolveMatrix(const TPZSolveMatrix &copy) {
        fTotalStrain = copy.fTotalStrain;
        fPlasticStrain = copy.fPlasticStrain;
        fSolution = copy.fSolution;
        fCmesh = copy.fCmesh;
        fNpts = copy.fNpts;
        fNphis = copy.fNphis;
        fStorage = copy.fStorage;
        fColSizes = copy.fColSizes;
        fRowSizes = copy.fRowSizes;
        fMatrixPosition = copy.fMatrixPosition;
        fColFirstIndex = copy.fColFirstIndex;
        fRowFirstIndex = copy.fRowFirstIndex;
        fElemColor = copy.fElemColor;
        fIndexes = copy.fIndexes;
        fIndexesColor = copy.fIndexesColor;
        fWeight = copy.fWeight;
        fMaterialData = copy.fMaterialData;

#ifdef __CUDACC__
        d_fStorage = copy.d_fStorage;
        d_fColSizes = copy.d_fColSizes;
        d_fRowSizes = copy.d_fRowSizes;
        d_fMatrixPosition = copy.d_fMatrixPosition;
        d_fColFirstIndex = copy.d_fColFirstIndex;
        d_fRowFirstIndex = copy.d_fRowFirstIndex;
        d_fElemColor = copy.d_fElemColor;
        d_fIndexes = copy.d_fIndexes;
        d_fIndexesColor = copy.d_fIndexesColor;

        d_GlobalSolution = copy.d_GlobalSolution;
        d_ExpandSolution = copy.d_ExpandSolution;
        d_Result = copy.d_Result;
        d_Weight = copy.d_Weight;
        d_Sigma = copy.d_Sigma;
        d_NodalForces = copy.d_NodalForces;
        d_GlobalForces = copy.d_GlobalForces;
#endif
    }

    TPZSolveMatrix &operator=(const TPZSolveMatrix &copy) {
        fTotalStrain = copy.fTotalStrain;
        fPlasticStrain = copy.fPlasticStrain;
        fSolution = copy.fSolution;
        fCmesh = copy.fCmesh;
        fStorage = copy.fStorage;
        fColSizes = copy.fColSizes;
        fRowSizes = copy.fRowSizes;
        fMatrixPosition = copy.fMatrixPosition;
        fColFirstIndex = copy.fColFirstIndex;
        fRowFirstIndex = copy.fRowFirstIndex;
        fElemColor = copy.fElemColor;
        fIndexes = copy.fIndexes;
        fIndexesColor = copy.fIndexesColor;
        fWeight = copy.fWeight;
        fMaterialData = copy.fMaterialData;

#ifdef __CUDACC__
        d_fStorage = copy.d_fStorage;
        d_fColSizes = copy.d_fColSizes;
        d_fRowSizes = copy.d_fRowSizes;
        d_fMatrixPosition = copy.d_fMatrixPosition;
        d_fColFirstIndex = copy.d_fColFirstIndex;
        d_fRowFirstIndex = copy.d_fRowFirstIndex;
        d_fElemColor = copy.d_fElemColor;
        d_fIndexes = copy.d_fIndexes;
        d_fIndexesColor = copy.d_fIndexesColor;
        d_Weight = copy.d_Weight;

        d_GlobalSolution = copy.d_GlobalSolution;
        d_ExpandSolution = copy.d_ExpandSolution;
        d_Result = copy.d_Result;
        d_Sigma = copy.d_Sigma;
        d_NodalForces = copy.d_NodalForces;
        d_GlobalForces = copy.d_GlobalForces;
#endif
        return *this;
    }

    void SetRowandColSizes(TPZVec<int> rowsize, TPZVec<int> colsize) {
        int64_t nelem = rowsize.size();

        fRowSizes.resize(nelem);
        fColSizes.resize(nelem);
        fMatrixPosition.resize(nelem + 1);
        fRowFirstIndex.resize(nelem + 1);
        fColFirstIndex.resize(nelem + 1);

        fRowSizes = rowsize;
        fColSizes = colsize;

        fMatrixPosition[0] = 0;
        fRowFirstIndex[0] = 0;
        fColFirstIndex[0] = 0;

        for (int64_t iel = 0; iel < nelem; ++iel) {
            fMatrixPosition[iel + 1] = fMatrixPosition[iel] + fRowSizes[iel] * fColSizes[iel];
            fRowFirstIndex[iel + 1] = fRowFirstIndex[iel] + fRowSizes[iel];
            fColFirstIndex[iel + 1] = fColFirstIndex[iel] + fColSizes[iel];
        }
        fStorage.resize(fMatrixPosition[nelem]);
        fElemColor.resize(nelem);
        fElemColor.Fill(-1);
    }

    void SetElementMatrix(int iel, TPZFMatrix<REAL> &elmat) {
        int64_t rows = fRowSizes[iel];
        int64_t cols = fColSizes[iel];
        int64_t pos = fMatrixPosition[iel];
        int64_t nelem = fRowSizes.size();

        TPZFMatrix<REAL> elmatloc(rows, cols, &fStorage[pos], rows*cols);
        elmatloc = elmat;
    }

    void SetIndexes(TPZVec<MKL_INT> indexes) {
        int64_t indsize = indexes.size();
        fIndexes.resize(indsize);
        fIndexes = indexes;
        fIndexesColor.resize(indsize);
    }

    void SetNumberofIntPoints(int64_t npts) {
        fNpts = npts;
    }

    int64_t NumberofIntPoints() {
        return fNpts;
    }

    void SetNumberofPhis(int64_t nphis) {
        fNphis = nphis;
    }

    void SetCompMesh(TPZCompMesh *cmesh) {
        fCmesh = cmesh;
    }

    void SetWeightVector (TPZStack<REAL> wvec) {
        fWeight = wvec;
    }

    void SetMaterialData (TElastoPlasticData materialdata) {
        fMaterialData = materialdata;
    }

    void SetSolution (TPZFMatrix<REAL> sol) {
        fSolution = sol;
    }

    void SetTotalStrain (TPZFMatrix<REAL> totalstrain) {
        fTotalStrain = totalstrain;
    }

    void SetPlasticStrain (TPZFMatrix<REAL> plasticstrain) {
        fPlasticStrain = plasticstrain;
    }

    void SetDataStructure();

    void GatherSolution(TPZFMatrix<REAL> &global_solution, TPZFMatrix<REAL> &gather_solution);

    void DeltaStrain(TPZFMatrix<REAL> &global_solution, TPZFMatrix<REAL> &deltastrain);

    void TotalStrain(TPZFMatrix<REAL> &delta_strain, TPZFMatrix<REAL> &total_strain);
    void ElasticStrain(TPZFMatrix<REAL> &total_strain, TPZFMatrix<REAL> &plastic_strain, TPZFMatrix<REAL> &elastic_strain);
    void PlasticStrain(TPZFMatrix<REAL> &total_strain, TPZFMatrix<REAL> &elastic_strain, TPZFMatrix<REAL> &plastic_strain);

    void ComputeStress(TPZFMatrix<REAL> &elastic_strain, TPZFMatrix<REAL> &sigma);

    void SpectralDecomposition(TPZFMatrix<REAL> &sigma_trial, TPZFMatrix<REAL> &eigenvalues, TPZFMatrix<REAL> &eigenvectors);
    void Normalize(double *sigma, double &maxel);
    void Interval(double *sigma, double *interval);
    void NewtonIterations(double *interval, double *sigma, double *eigenvalues, double &maxel);
    void Eigenvectors(double *sigma, double *eigenvalue, double *eigenvector, double &maxel);
    void Multiplicity1(double *sigma, double eigenvalue, double *eigenvector);
    void Multiplicity2(double *sigma, double eigenvalue, double *eigenvector1, double *eigenvector2);

    void ProjectSigma(TPZFMatrix<REAL> &eigenvalues, TPZFMatrix<REAL> &sigma_projected);
    bool PhiPlane(double *eigenvalues, double *sigma_projected);
    bool ReturnMappingMainPlane(double *eigenvalues, double *sigma_projected, double &m_hardening);
    bool ReturnMappingRightEdge(double *eigenvalues, double *sigma_projected, double &m_hardening);
    bool ReturnMappingLeftEdge(double *eigenvalues, double *sigma_projected, double &m_hardening);
    void ReturnMappingApex(double *eigenvalues, double *sigma_projected, double &m_hardening);

    void StressCompleteTensor(TPZFMatrix<REAL> &sigma_projected, TPZFMatrix<REAL> &eigenvectors, TPZFMatrix<REAL> &sigma);

    void ComputeStrain(TPZFMatrix<REAL> &sigma, TPZFMatrix<REAL> &elastic_strain);

    void NodalForces(TPZFMatrix<REAL> &sigma, TPZFMatrix<REAL> &nodal_forces);

    void ColoredAssemble(TPZFMatrix<REAL> &nodal_forces_vec, TPZFMatrix<REAL> &nodal_forces_global);

    TPZFMatrix<REAL> AssembleResidual();

    void ColoringElements() const;

protected:

///total strain
    TPZFMatrix<REAL> fTotalStrain;

///plastic strain
    TPZFMatrix<REAL> fPlasticStrain;

/// solution
    TPZFMatrix<REAL> fSolution;

/// computational mesh
    TPZCompMesh *fCmesh;

///total number of int points
    int64_t fNpts;

///total number of phis
    int64_t fNphis;

/// vector containing the matrix coefficients
    TPZVec<REAL> fStorage;

/// number of rows of each block matrix
    TPZVec<int> fRowSizes;

/// number of columns of each block matrix
    TPZVec<int> fColSizes;

/// indexes vector in x and y direction
    TPZVec<MKL_INT> fIndexes;

/// indexes vector in x and y direction by color
    TPZVec<MKL_INT> fIndexesColor;

/// color indexes of each element
    TPZVec<int64_t> fElemColor;

/// position of the matrix of the elements
    TPZVec<int> fMatrixPosition;

/// position of the result vector
    TPZVec<int> fRowFirstIndex;

/// position in the fIndex vector of each element
    TPZVec<int> fColFirstIndex;

/// Weight Vector
    TPZStack<REAL> fWeight;

/// material data
    TElastoPlasticData fMaterialData;

/// Parameters stored on device
#ifdef __CUDACC__
    REAL *d_fStorage;
    int *d_fColSizes;
    int *d_fRowSizes;
    int *d_fMatrixPosition;
    int *d_fColFirstIndex;
    int *d_fRowFirstIndex;
    int *d_fElemColor;
    int *d_fIndexes;
    int *d_fIndexesColor;

    REAL *d_GlobalSolution;
    REAL *d_ExpandSolution;
    REAL *d_Result;
    REAL *d_Weight;
    REAL *d_Sigma;
    REAL *d_NodalForces;
    REAL *d_GlobalForces;

    //library handles
    cusparseHandle_t handle_cusparse;
    cublasHandle_t handle_cublas;
#endif

};

#endif /* TPZSolveMatrix_h */
