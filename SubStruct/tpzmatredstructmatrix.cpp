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

template< class TStructMatrix, class TSparseMatrix>
TPZMatRedStructMatrix<TStructMatrix,TSparseMatrix>::TPZMatRedStructMatrix() : TPZStructMatrix()
{
}

template< class TStructMatrix, class TSparseMatrix>
void TPZMatRedStructMatrix<TStructMatrix,TSparseMatrix>::SetMesh(TPZCompMesh *cmesh) {
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


template< class TStructMatrix, class TSparseMatrix>
TPZMatRedStructMatrix<TStructMatrix,TSparseMatrix>::TPZMatRedStructMatrix(TPZSubCompMesh *mesh) : TPZStructMatrix(mesh)
{
	fInternalEqs = mesh->NumInternalEquations();
}

template< class TStructMatrix, class TSparseMatrix>
TPZMatRedStructMatrix<TStructMatrix,TSparseMatrix>::~TPZMatRedStructMatrix()
{
}

template< class TStructMatrix, class TSparseMatrix>
TPZMatRedStructMatrix<TStructMatrix,TSparseMatrix>::TPZMatRedStructMatrix(const TPZMatRedStructMatrix<TStructMatrix,TSparseMatrix> &copy) : TPZStructMatrix(copy)
{
	fInternalEqs = copy.fInternalEqs;
}

template< class TStructMatrix, class TSparseMatrix>
TPZStructMatrix *TPZMatRedStructMatrix<TStructMatrix,TSparseMatrix>::Clone()
{
	return new TPZMatRedStructMatrix(*this);
}

template< class TStructMatrix, class TSparseMatrix>
TPZMatrix<STATE> *TPZMatRedStructMatrix<TStructMatrix,TSparseMatrix>::Create()
{
	int neq = fMesh->NEquations();
	TStructMatrix strmat(fMesh);
	strmat.SetEquationRange(0,fInternalEqs);
	TPZMatrix<STATE> *InternalStiff = strmat.Create();
	TPZMatRed<STATE, TSparseMatrix> *matred = new TPZMatRed<STATE, TSparseMatrix>(neq,fInternalEqs);
	matred->SetK00(InternalStiff);
	return matred;
}

template< class TStructMatrix, class TSparseMatrix>
int TPZMatRedStructMatrix<TStructMatrix,TSparseMatrix>::ClassId() const{
    return Hash("TPZMatRedStructMatrix") ^ TPZStructMatrix::ClassId() << 1 ^
        TSparseMatrix().ClassId() << 2;
}

template< class TStructMatrix, class TSparseMatrix>
void TPZMatRedStructMatrix<TStructMatrix,TSparseMatrix>::Read(TPZStream& buf, void* context){
    TPZStructMatrix::Read(buf,context);
    buf.Read(&fInternalEqs);
}

template< class TStructMatrix, class TSparseMatrix>
void TPZMatRedStructMatrix<TStructMatrix,TSparseMatrix>::Write(TPZStream& buf, int withclassid) const{
    TPZStructMatrix::Write(buf,withclassid);
    buf.Write(fInternalEqs);
}


template class TPZMatRedStructMatrix<TPZSkylineStructMatrix<STATE>,TPZVerySparseMatrix<STATE> >;
template class TPZMatRedStructMatrix<TPZSkylineStructMatrix<STATE>,TPZFMatrix<STATE> >;