/**
 * @file
 * @brief Contains the TPZSolveMatrix class which implements a solution based on a matrix procedure.
 */

#ifndef TPZIntPointsFEM_h
#define TPZIntPointsFEM_h

//#include <StrMatrix/TPZSSpStructMatrix.h>
#include <TPZSSpStructMatrix.h>
#include "pzsysmp.h"
#include "tpzverysparsematrix.h"
#include "TPZIrregularBlocksMatrix.h"
#include "TPZNumericalIntegrator.h"
#include <map>


class TPZElastoPlasticIntPointsStructMatrix : public TPZSymetricSpStructMatrix {
    
public:
    
    /** @brief Default constructor */
    TPZElastoPlasticIntPointsStructMatrix();

    /** @brief Creates the object based on a Compmesh
     * @param Compmesh : Computational mesh */
    TPZElastoPlasticIntPointsStructMatrix(TPZCompMesh *cmesh);

    /** @brief Default destructor */
    ~TPZElastoPlasticIntPointsStructMatrix();

    /** @brief Clone */
    TPZStructMatrix *Clone();
    
    TPZMatrix<STATE> * Create();

    // need help
    TPZMatrix<STATE> *CreateAssemble(TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);

    void SetUpDataStructure();

    void Assemble(TPZMatrix<STATE> & mat, TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
    
    void Assemble(TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
    
    bool isBuilt() {
        if(fIntegrator.IrregularBlocksMatrix().Rows() != 0) return true;
        else return false;
    }

#ifdef USING_CUDA
    void TransferDataToGPU();
#endif

private:
    
    int StressRateVectorSize();
    
    void ComputeDomainElementIndexes(TPZVec<int> &element_indexes);
    
    void ClassifyMaterialsByDimension();
    
    void AssembleBoundaryData();

    void SetUpIrregularBlocksData(TPZVec<int> &element_indexes, TPZIrregularBlocksMatrix::IrregularBlocks &blocksData);

    void SetUpIndexes(TPZVec<int> &element_indexes, TPZVec<int> & dof_indexes);

    void ColoredIndexes(TPZVec<int> &element_indexes, TPZVec<int> &indexes, TPZVec<int> &coloredindexes, int &ncolor);
    
    void FillLIndexes();
    
    int fDimension;
    
    int64_t me(TPZVec<int> &IA, TPZVec<int> &JA, int64_t & i_dest, int64_t & j_dest);
    
    TPZNumericalIntegrator fIntegrator;

    TPZVerySparseMatrix<STATE> fSparseMatrixLinear; //-> BC data

    TPZFMatrix<STATE> fRhsLinear; //-> BC data
    
    std::set<int> fBCMaterialIds;
    
    TPZVec<int> m_IA_to_sequence;
    
    TPZVec<int> m_JA_to_sequence;
    
    TPZVec<int> m_color_l_sequence;
    
    TPZVec<int> m_first_color_l_index;
    
    std::vector<int64_t> m_el_color_indexes;
    
    std::vector<int64_t> m_first_color_index;

    #ifdef USING_CUDA
    TPZCudaCalls fCudaCalls;

    TPZVecGPU<int> d_color_l_sequence;

    TPZVecGPU<int> d_IA_to_sequence;
    
    TPZVecGPU<int> d_JA_to_sequence;
    
    TPZVecGPU<int64_t> d_el_color_indexes;

    TPZVecGPU<REAL> d_RhsLinear;
    #endif
    
};

#endif /* TPZIntPointsFEM_h */
