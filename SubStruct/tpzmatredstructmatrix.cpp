/**
 * @file
 * @brief Contains the implementation of the TPZMatRedStructMatrix methods. 
 * @since 22/04/2009
 */

#include "tpzmatredstructmatrix.h"
#include "pzskylstrmatrix.h"
#include "tpzverysparsematrix.h"
#include "pzsubcmesh.h"
#include "pzmatred.h"

template< class TStructMatrix, class TSparseMatrix, class TVar, class TPar>
TPZMatRedStructMatrix<TStructMatrix,TSparseMatrix,TVar,TPar>::TPZMatRedStructMatrix() : TPZStructMatrixT<TVar>()
{
}

template< class TStructMatrix, class TSparseMatrix, class TVar, class TPar>
void TPZMatRedStructMatrix<TStructMatrix,TSparseMatrix,TVar,TPar>::SetMesh(TPZCompMesh *cmesh) {
    TPZStructMatrix::SetMesh(cmesh);
    if (cmesh){
        TPZSubCompMesh *subCMesh = dynamic_cast<TPZSubCompMesh *>(cmesh);
        if (!subCMesh){
            DebugStop();
        }
	fInternalEqs = subCMesh->NumInternalEquations();
    } else {
        DebugStop();
    }
}


template< class TStructMatrix, class TSparseMatrix, class TVar, class TPar>
TPZMatRedStructMatrix<TStructMatrix,TSparseMatrix,TVar,TPar>::TPZMatRedStructMatrix(TPZSubCompMesh *mesh) : TPZStructMatrixT<TVar>(mesh)
{
	fInternalEqs = mesh->NumInternalEquations();
}

template< class TStructMatrix, class TSparseMatrix, class TVar, class TPar>
TPZMatRedStructMatrix<TStructMatrix,TSparseMatrix,TVar,TPar>::~TPZMatRedStructMatrix()
{
}

template< class TStructMatrix, class TSparseMatrix, class TVar, class TPar>
TPZMatRedStructMatrix<TStructMatrix,TSparseMatrix,TVar,TPar>::TPZMatRedStructMatrix(const TPZMatRedStructMatrix<TStructMatrix,TSparseMatrix,TVar,TPar> &copy) : TPZStructMatrixT<TVar>(copy)
{
	fInternalEqs = copy.fInternalEqs;
}

template< class TStructMatrix, class TSparseMatrix, class TVar, class TPar>
TPZStructMatrix *TPZMatRedStructMatrix<TStructMatrix,TSparseMatrix,TVar,TPar>::Clone()
{
	return new TPZMatRedStructMatrix(*this);
}

template< class TStructMatrix, class TSparseMatrix, class TVar, class TPar>
TPZMatrix<STATE> *TPZMatRedStructMatrix<TStructMatrix,TSparseMatrix,TVar,TPar>::Create()
{
	int neq = this->fMesh->NEquations();
	TStructMatrix strmat(this->fMesh);
	strmat.SetEquationRange(0,fInternalEqs);
	TPZMatrix<STATE> *InternalStiff = strmat.Create();
	TPZMatRed<STATE, TSparseMatrix> *matred = new TPZMatRed<STATE, TSparseMatrix>(neq,fInternalEqs);
	matred->SetK00(InternalStiff);
	return matred;
}

template< class TStructMatrix, class TSparseMatrix, class TVar, class TPar>
int TPZMatRedStructMatrix<TStructMatrix,TSparseMatrix,TVar,TPar>::ClassId() const{
    return Hash("TPZMatRedStructMatrix") ^
        TPZStructMatrix::ClassId() << 1 ^
        TSparseMatrix().ClassId() << 2 ^
        TPar::ClassId() << 3;
}

template< class TStructMatrix, class TSparseMatrix, class TVar, class TPar>
void TPZMatRedStructMatrix<TStructMatrix,TSparseMatrix,TVar,TPar>::Read(TPZStream& buf, void* context){
    TPZStructMatrix::Read(buf,context);
    buf.Read(&fInternalEqs);
}

template< class TStructMatrix, class TSparseMatrix, class TVar, class TPar>
void TPZMatRedStructMatrix<TStructMatrix,TSparseMatrix,TVar,TPar>::Write(TPZStream& buf, int withclassid) const{
    TPZStructMatrix::Write(buf,withclassid);
    buf.Write(fInternalEqs);
}


template class TPZMatRedStructMatrix<TPZSkylineStructMatrix<STATE>,TPZVerySparseMatrix<STATE>,STATE,TPZStructMatrixOR<STATE> >;
template class TPZMatRedStructMatrix<TPZSkylineStructMatrix<STATE>,TPZFMatrix<STATE>,STATE,TPZStructMatrixOR<STATE> >;