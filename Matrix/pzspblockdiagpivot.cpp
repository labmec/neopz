/**
 * @file
 * @brief Contains the implementation of the TPZSpBlockDiagPivot methods.
 */
//
// C++ Implementation: %{MODULE}
//
// Description:
//
//
// Author: %{AUTHOR} <%{EMAIL}>, (C) %{YEAR}
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "pzspblockdiagpivot.h"
#include "pzfmatrix.h" 
using namespace std;

TPZSpBlockDiagPivot::TPZSpBlockDiagPivot()
: TPZSparseBlockDiagonal(), fPivotIndices()
{}


TPZSpBlockDiagPivot::~TPZSpBlockDiagPivot()
{}

int TPZSpBlockDiagPivot::Decompose_LU(){
	if (  fDecomposed && fDecomposed == ELU) {
		return ELU;
	} else if(fDecomposed) {
		PZError << __PRETTY_FUNCTION__ << endl;
		PZError << "TPZSpBlockDiagPivot::Decompose_LU() - Matrix is already decomposed with other scheme" << endl;
	}
	
	this->fPivotIndices.Resize( fBlock.NElements() );
#ifdef DEBUG  
	this->fPivotIndices.Fill( -1 );
#endif
	TPZManVector<int> pivot;
	int pivotindex = 0;
	const int nb = fBlockSize.NElements();
	for(int b = 0; b < nb; b++) {
		const int pos = fBlockPos[b];
		const int bsize = fBlockSize[b];
		if(!bsize) continue;
		TPZFMatrix temp(bsize,bsize,&fStorage[pos],bsize*bsize);
		temp.Decompose_LU(pivot);
		for(int id = 0; id < bsize; id++){
			this->fPivotIndices[pivotindex + id] = pivot[id];
		}
		pivotindex += bsize;
	}
	fDecomposed = ELUPivot;
#ifdef DEBUG
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

int TPZSpBlockDiagPivot::Substitution( TPZFMatrix * B ) const{
	TPZFNMatrix<1000> BG(fBlock.NElements(),B->Cols());
	Gather(*B,BG,1);
	int result = this->Substitution2(&BG);
	B->Zero();
	Scatter(BG,*B,1);
	return result;
}

int TPZSpBlockDiagPivot::Substitution2( TPZFMatrix *B) const
{
	if(fDecomposed != ELUPivot) {
		PZError << __PRETTY_FUNCTION__ << "TPZSpBlockDiagPivot::Decompose_LU is decomposed with other scheme than ELUPivot." << endl;
	}
	
	TPZManVector<int, 1000> pivot;
	int b,eq=0;
	const int nb = fBlockSize.NElements();
	int c, nc = B->Cols();
	for(c=0; c<nc; c++) {
		eq = 0;
		int pivotindex = 0;
		for(b=0;b<nb; b++) {
			//      const int pos = fBlockPos[b];
			const int bsize = fBlockSize[b];
			if(!bsize) continue;
			//      TPZFMatrix temp(bsize,bsize,&fStorage[pos],bsize*bsize);
			//      temp.SetIsDecomposed(ELUPivot);
			TPZFMatrix BTemp(bsize,1,&(B->operator()(eq,c)),bsize);
			pivot.Resize(bsize);
			//      memcpy(&pivot[0],&fPivotIndices[pivotindex],bsize*sizeof(int));
			for(int id = 0; id < bsize; id++){
				pivot[id] = this->fPivotIndices[pivotindex+id];
			}
			pivotindex += bsize;
			TPZFMatrix::Substitution(fStorage,bsize,&BTemp);
			//      temp.Substitution(&BTemp, pivot);
			eq+= bsize;
		}
	}
	return 1;
}


