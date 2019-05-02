/**
 * @file
 * @brief Contains the TPZSolveMatrix class which implements a solution based on a matrix procedure.
 */

#ifndef TPZIntPointsFEM_h
#define TPZIntPointsFEM_h
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzinterpolationspace.h"
#include "pzcmesh.h"
#include "TPZPlasticStepPV.h"
#include "TPZYCMohrCoulombPV.h"
#include "TPZElastoPlasticMem.h"
#include "TPZMatElastoPlastic2D.h"

#ifdef USING_MKL
#include "mkl.h"
#endif


#ifdef __CUDACC__
//#include "TPZVecGPU.h"
#include <cublas_v2.h>
#include <cusparse.h>
#include <cuda.h>
#endif

class TPZIntPointsFEM {

public:

    TPZIntPointsFEM();

    TPZIntPointsFEM(TPZCompMesh *cmesh, int materialid);

    ~TPZIntPointsFEM();

    TPZIntPointsFEM(const TPZIntPointsFEM &copy);

    TPZIntPointsFEM &operator=(const TPZIntPointsFEM &copy);

    void SetRowandColSizes(TPZVec<int> rowsize, TPZVec<int> colsize) {
        int64_t nelem = rowsize.size();

        fRowSizes.resize(nelem);
        fColSizes.resize(nelem);
        fMatrixPosition.resize(nelem + 1);
        fRowFirstIndex.resize(nelem + 1);
        fColFirstIndex.resize(nelem + 1);
        fMatrixPosition[0] = 0;
        fRowFirstIndex[0] = 0;
        fColFirstIndex[0] = 0;

        for (int64_t iel = 0; iel < nelem; ++iel) {
        	fRowSizes[iel] = rowsize[iel];
        	fColSizes[iel] = colsize[iel];
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

        TPZFMatrix<REAL> elmatloc(rows, cols, &fStorage[pos], rows*cols);
        elmatloc = elmat;
    }

    void SetIndexes(TPZVec<int> indexes) {
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

    void SetWeightVector (TPZVec<REAL> &wvec) {
        fWeight = wvec;
    }

    void SetMaterialId (int materialid) {
        TPZMaterial *material = fCmesh->FindMaterial(materialid);
        fMaterial = dynamic_cast<TPZMatElastoPlastic2D<TPZPlasticStepPV<TPZYCMohrCoulombPV,TPZElasticResponse> , TPZElastoPlasticMem> *>(material);
    }

    void LoadSolution (TPZFMatrix<REAL> & sol) {
        fSolution = sol;
    }

     TPZFMatrix<REAL> & Rhs() {
        return fRhs;
    }

    void SetMeshDimension(int dim) {
        fDim = dim;
    }

    void SetDataStructure();

    void GatherSolution(TPZFMatrix<REAL> &global_solution, TPZFMatrix<REAL> &gather_solution);
    void DeltaStrain(TPZFMatrix<REAL> &gather_solution, TPZFMatrix<REAL> &delta_strain);

    void ElasticStrain(TPZFMatrix<REAL> &delta_strain, TPZFMatrix<REAL> &plastic_strain, TPZFMatrix<REAL> &elastic_strain);
    void PlasticStrain(TPZFMatrix<REAL> &delta_strain, TPZFMatrix<REAL> &elastic_strain, TPZFMatrix<REAL> &plastic_strain);

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

    void ColoredAssemble(TPZFMatrix<REAL> &nodal_forces, TPZFMatrix<REAL> &residual);

    void AssembleResidual();

    void ColoringElements() const;

    void AssembleRhsBoundary();

#ifdef __CUDACC__
    void TransferDataStructure();

    void GatherSolutionGPU(REAL *gather_solution);
    void DeltaStrainGPU(REAL *gather_solution, REAL *delta_strain);
    void ElasticStrainGPU(REAL *delta_strain, REAL *plastic_strain, REAL *elastic_strain);
    void ComputeStressGPU(REAL *elastic_strain, REAL *sigma);
    void SpectralDecompositionGPU(REAL *sigma_trial, REAL *eigenvalues, REAL *eigenvectors);
    void ProjectSigmaGPU(REAL *eigenvalues, REAL *sigma_projected);
    void StressCompleteTensorGPU(REAL *sigma_projected, REAL *eigenvectors, REAL *sigma);
    void NodalForcesGPU(REAL *sigma, REAL *nodal_forces);
    void ColoredAssembleGPU(REAL *nodal_forces, REAL *residual);
    void ComputeStrainGPU(REAL *sigma, REAL *elastic_strain);
    void PlasticStrainGPU(REAL *delta_strain, REAL *elastic_strain, REAL *plastic_strain);

#endif

protected:
    int fDim;
    TPZStack<int64_t> fBoundaryElements;
    TPZCompMesh *fCmesh;
    int64_t fNpts;
    int64_t fNphis;
    TPZVec<int64_t> fElemColor;
    TPZMatElastoPlastic2D<TPZPlasticStepPV<TPZYCMohrCoulombPV,TPZElasticResponse>, TPZElastoPlasticMem> *fMaterial;

    TPZFMatrix<REAL> fRhs;
    TPZFMatrix<REAL> fRhsBoundary;
	TPZFMatrix<REAL> fSolution;
    TPZFMatrix<REAL> fPlasticStrain;
	TPZVec<REAL> fStorage;
	TPZVec<int> fRowSizes;
	TPZVec<int> fColSizes;
	TPZVec<int> fMatrixPosition;
	TPZVec<int> fRowFirstIndex;
	TPZVec<int> fColFirstIndex;
	TPZVec<int> fIndexes;
	TPZVec<int> fIndexesColor;
	TPZVec<REAL> fWeight;

#ifdef __CUDACC__
	cusparseHandle_t handle_cusparse;
	cublasHandle_t handle_cublas;
#endif

//	 @omar:: Cara Natalia AO DEFINIR dSolution DA ERRADO COMENTA dSolution QUALQUER COISA NO wapp

    REAL *dRhs;
    REAL *dRhsBoundary;
    REAL *dSolution;
    REAL *dPlasticStrain;
    REAL *dStorage;
    int *dRowSizes;
    int *dColSizes;
    int *dMatrixPosition;
    int *dRowFirstIndex;
    int *dColFirstIndex;
    int *dIndexes;
    int *dIndexesColor;
    REAL *dWeight;
//#endif


};

#endif /* TPZIntPointsFEM_h */
