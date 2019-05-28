/**
 * @file
 * @brief Contains the TPZSolveMatrix class which implements a solution based on a matrix procedure.
 */

#ifndef TPZIntPointsFEM_h
#define TPZIntPointsFEM_h
#include "TPZIrregularBlockMatrix.h"
#include "TPZIntPointsData.h"
#include "pzstrmatrix.h"

class TPZIntPointsStructMatrix : public TPZStructMatrix {
public:
    /** @brief Default constructor */
    TPZIntPointsStructMatrix();

    /** @brief Creates the object based on a Compmesh
     * @param Compmesh : Computational mesh
     */
    TPZIntPointsStructMatrix(TPZCompMesh *cmesh);

    TPZIntPointsStructMatrix(TPZAutoPointer<TPZCompMesh> cmesh);

    /** @brief Default destructor */
    ~TPZIntPointsStructMatrix();

    /** @brief Clone */
    TPZStructMatrix *Clone();

    /** @brief Creates a TPZIntPointsStructMatrix object with copy constructor
     * @param copy : original TPZIntPointsStructMatrix object
     */
    TPZIntPointsStructMatrix(const TPZIntPointsStructMatrix &copy);

    /** @brief operator= */
    TPZIntPointsStructMatrix &operator=(const TPZIntPointsStructMatrix &copy);

    /** @brief Defines with elements must be assembled */
    void ElementsToAssemble();

    /** @brief Defines matrices rows and columns sizes, first row and column indexes and CSR parameters */
    void BlocksInfo(int imat);

    /** @brief Fill matrices values */
    void FillBlocks(int imat);

    /** @brief Defines integration points information */
    void IntPointsInfo(int imat);

    /** @brief Assemble the load vector */
    void Assemble(TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);

    /** @brief Performs elements coloring */
    void ColoringElements(int imat);

    /** @brief Assemble the load vector for boundary elements */
    void AssembleRhsBoundary();

    /** @brief Initialize fBlockMatrix and fIntPointsData */
    void Initialize();

    /** @brief Access methods */
    TPZIrregularBlockMatrix BlockMatrix(int imat) {
        return fBlockMatrix[imat];
    }

    TPZIntPointsData IntPointsData(int imat) {
        return fIntPointsData[imat];
    }
    TPZFMatrix<REAL> Rhs() {
        return fRhs;
    }

protected:
    /** @brief Load vector */
    TPZFMatrix<REAL> fRhs;

    /** @brief Load vector for boundary elements */
    TPZFMatrix<REAL> fRhsBoundary;

    /** @brief Vector of irregular blocks matrix (each position of the vector represents one material) */
    TPZVec<TPZIrregularBlockMatrix> fBlockMatrix;

    /** @brief Vector of integration points info (each position of the vector represents one material) */
    TPZVec<TPZIntPointsData> fIntPointsData;

    /** @brief Vector of indexes of the domain elements (each position of the vector represents one material) */
    TPZVec<TPZStack<int64_t>> fElemIndexes;
};

#endif /* TPZIntPointsFEM_h */
