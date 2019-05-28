/**
 * @file
 * @brief Contains the TPZSolveMatrix class which implements a solution based on a matrix procedure.
 */

#ifndef TPZIntPointsFEM_h
#define TPZIntPointsFEM_h

#include <StrMatrix/TPZSSpStructMatrix.h>
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

    // need help
    TPZMatrix<STATE> *CreateAssemble(TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);

    void AssembleRhsBoundary(TPZFMatrix<REAL> &rhsboundary);

    void SetUpDataStructure();

    void CalcResidual(TPZFMatrix<REAL> &rhs);

    bool isBuilt() {
        if(fCoefToGradSol.IrregularBlocksMatrix().Rows() != 0) return true;
        else return false;
    }

private:
    TPZCoefToGradSol fCoefToGradSol;

    TPZMyLambdaExpression fLambdaExp;

    TPZStructMatrix fStructMatrix;
};

#endif /* TPZIntPointsFEM_h */
