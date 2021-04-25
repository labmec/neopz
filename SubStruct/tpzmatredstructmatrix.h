/**
 * @file
 * @brief Contains the TPZMatRedStructMatrix class. 
 */

#ifndef TPZMATREDSTRUCTMATRIX
#define TPZMATREDSTRUCTMATRIX

#include "TPZStructMatrixT.h"
#include "pzstrmatrixor.h"
class TPZSubCompMesh;


/**
 * @ingroup substructure
 * @brief .. . \ref substructure "Sub Structure"
 */
template<
    class TStructMatrix, class TSparseMatrix,
    class TVar=STATE, class TPar=TPZStructMatrixOR<TVar>
    >
class TPZMatRedStructMatrix : public TPZStructMatrixT<TVar>,
                              public TPar
{
public:
	/** @brief Constructor */
	TPZMatRedStructMatrix(TPZSubCompMesh *mesh);
	/** @brief Destructor */
	virtual ~TPZMatRedStructMatrix();
	/** @brief Copy constructor */
	TPZMatRedStructMatrix(const TPZMatRedStructMatrix &copy);
	
	TPZStructMatrix *Clone() override;
	
	TPZMatrix<STATE> *Create() override;

	//@{
    //!Read and Write methods
    int ClassId() const override;

    void Read(TPZStream& buf, void* context) override;

    void Write(TPZStream& buf, int withclassid) const override;
    //@}
        
protected :
	void SetMesh(TPZCompMesh *cmesh);
        
private:
	TPZMatRedStructMatrix();
    
	friend TPZPersistenceManager;
	
	int fInternalEqs;
	
};

#endif
