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
#include "TPZMyLambdaExpression.h"
#include "TPZCoefToGradSol.h"


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

    void AssembleBoundaryData(std::set<int> &boundary_matids);

    void SetUpDataStructure();

    void Assemble(TPZMatrix<STATE> & mat, TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
    
    void Assemble(TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
    
    bool isBuilt() {
        if(fCoefToGradSol.IrregularBlocksMatrix().Rows() != 0) return true;
        else return false;
    }

    void Dep(TPZVec<REAL> &depxx, TPZVec<REAL> &depyy, TPZVec<REAL> &depxy);

private:
    void GetDomainElements(TPZStack<int> &elindex_domain, std::set<int> &boundary_matids);

    void SetUpIrregularBlocksData(TPZStack<int> &elindex_domain, TPZIrregularBlocksMatrix::IrregularBlocks &blocksData);

    void SetUpIndexes(TPZStack<int> &elindex_domain, TPZVec<int> & dof_indexes);

    void ColoredIndexes(TPZStack<int> &elindex_domain, TPZVec<int> &indexes, TPZVec<int> &coloredindexes, int &ncolor);

    TPZCoefToGradSol fCoefToGradSol;

    TPZMyLambdaExpression fLambdaExp;

    TPZVerySparseMatrix<STATE> fSparseMatrixLinear; //-> BC data

    TPZFMatrix<STATE> fRhsLinear; //-> BC data
    
};

#endif /* TPZIntPointsFEM_h */
