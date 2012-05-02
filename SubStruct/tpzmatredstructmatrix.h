/**
 * @file
 * @brief Contains the TPZMatRedStructMatrix class. 
 */

#ifndef TPZMATREDSTRUCTMATRIX
#define TPZMATREDSTRUCTMATRIX

#include "pzstrmatrix.h"
class TPZSubCompMesh;
template <class TVar>
class TPZFMatrix;

/**
 * @ingroup substructure
 * @brief .. . \ref substructure "Sub Structure"
 */
template<class TStructMatrix, class TSparseMatrix>
class TPZMatRedStructMatrix : TPZStructMatrix
{
public:
	/** @brief Constructor */
	TPZMatRedStructMatrix(TPZSubCompMesh *mesh);
	/** @brief Destructor */
	virtual ~TPZMatRedStructMatrix();
	/** @brief Copy constructor */
	TPZMatRedStructMatrix(const TPZMatRedStructMatrix &copy);
	
	virtual TPZStructMatrix *Clone();
	
	virtual TPZMatrix<STATE> *Create();
	
private:
	
	int fInternalEqs;
	
};

#endif
