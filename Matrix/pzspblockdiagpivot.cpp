/**
 * @file
 * @brief Contains the implementation of the TPZSpBlockDiagPivot methods.
 */

#include "pzspblockdiagpivot.h"
#include "pzfmatrix.h" 
using namespace std;

template<class TVar>
TPZSpBlockDiagPivot<TVar>::TPZSpBlockDiagPivot()
: TPZSparseBlockDiagonal<TVar>(), fPivotIndices()
{}

template<class TVar>
TPZSpBlockDiagPivot<TVar>::~TPZSpBlockDiagPivot()
{}

template<class TVar>
int TPZSpBlockDiagPivot<TVar>::Decompose_LU(){
	if (  this->fDecomposed && this->fDecomposed == ELU) {
		return ELU;
	} else if(this->fDecomposed) {
		PZError << __PRETTY_FUNCTION__ << endl;
		PZError << "TPZSpBlockDiagPivot::Decompose_LU() - Matrix is already decomposed with other scheme" << endl;
	}
	
	this->fPivotIndices.Resize( this->fBlock.NElements() );
#ifdef PZDEBUG  
	this->fPivotIndices.Fill( -1 );
#endif
	TPZManVector<int> pivot;
	int pivotindex = 0;
	const int nb = this->fBlockSize.NElements();
	for(int b = 0; b < nb; b++) {
		const int pos = this->fBlockPos[b];
		const int bsize = this->fBlockSize[b];
		if(!bsize) continue;
		TPZFMatrix<TVar> temp(bsize,bsize,&this->fStorage[pos],bsize*bsize);
		temp.Decompose_LU(pivot);
		for(int id = 0; id < bsize; id++){
			this->fPivotIndices[pivotindex + id] = pivot[id];
		}
		pivotindex += bsize;
	}
	this->fDecomposed = ELUPivot;
#ifdef PZDEBUG
	{
		int n = this->fPivotIndices.NElements();
		for(int i = 0; i < n; i++){
			if ( this->fPivotIndices[i] == -1 ){
				PZError << __PRETTY_FUNCTION__ << endl;
				PZError << "TPZSpBlockDiagPivot::Decompose_LU() - fPivotIndices attribute has error" << endl;      
			}
		}      
	}
#endif  
	return 1;
}

template<class TVar>
int TPZSpBlockDiagPivot<TVar>::Substitution( TPZFMatrix<TVar> * B ) const{
	TPZFNMatrix<1000,TVar> BG(this->fBlock.NElements(),B->Cols());
	Gather(*B,BG,1);
	int result = this->Substitution2(&BG);
	B->Zero();
	Scatter(BG,*B,1);
	return result;
}

template<class TVar>
int TPZSpBlockDiagPivot<TVar>::Substitution2( TPZFMatrix<TVar> *B) const
{
	if(this->fDecomposed != ELUPivot) {
		PZError << __PRETTY_FUNCTION__ << "TPZSpBlockDiagPivot::Decompose_LU is decomposed with other scheme than ELUPivot." << endl;
	}
	
	TPZManVector<int, 1000> pivot;
	int b,eq=0;
	const int nb = this->fBlockSize.NElements();
	int c, nc = B->Cols();
	for(c=0; c<nc; c++) {
		eq = 0;
		int pivotindex = 0;
		for(b=0;b<nb; b++) {
			const int bsize = this->fBlockSize[b];
			if(!bsize) continue;
			TPZFMatrix<TVar> BTemp(bsize,1,&(B->operator()(eq,c)),bsize);
			pivot.Resize(bsize);
			for(int id = 0; id < bsize; id++){
				pivot[id] = this->fPivotIndices[pivotindex+id];
			}
			pivotindex += bsize;
			TPZFMatrix<TVar>::Substitution(this->fStorage,bsize,&BTemp);
			eq+= bsize;
		}
	}
	return 1;
}
