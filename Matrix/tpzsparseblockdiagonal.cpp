
/**
 * @file
 * @brief Contains the implementation of the TPZSparseBlockDiagonal methods.
 */

#include "tpzsparseblockdiagonal.h"
#include "pzfmatrix.h"

#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.StrMatrix");
#endif

using namespace std;

template<class TVar>
TPZSparseBlockDiagonal<TVar>::TPZSparseBlockDiagonal() : TPZRegisterClassId(&TPZSparseBlockDiagonal::ClassId)
{
}
template<class TVar>
TPZSparseBlockDiagonal<TVar>::TPZSparseBlockDiagonal(TPZVec<int64_t> &blockgraph, TPZVec<int64_t> &blockgraphindex,int64_t rows) : TPZRegisterClassId(&TPZSparseBlockDiagonal::ClassId), fBlock(blockgraph), fBlockIndex(blockgraphindex)
{
	int64_t numbl = blockgraphindex.NElements()-1;
	this->fBlockSize.Resize(numbl);
	int64_t ibl;
	for(ibl=0; ibl<numbl; ibl++)
	{
		this->fBlockSize[ibl] = blockgraphindex[ibl+1]-blockgraphindex[ibl];
	}
	this->Initialize(this->fBlockSize);
	this->fRow = rows;
	this->fCol = rows;
}

template<class TVar>
TPZSparseBlockDiagonal<TVar>::TPZSparseBlockDiagonal(TPZVec<int64_t> &blockgraph, TPZVec<int64_t> &blockgraphindex,int64_t rows, int color, TPZVec<int> &colors) : TPZRegisterClassId(&TPZSparseBlockDiagonal::ClassId)
{
#ifdef PZ_LOG
	LOGPZ_DEBUG(logger, "Constructor of TPZSparseBlockDiagonal");
#endif
	int64_t numbl = blockgraphindex.NElements()-1;
	this->fBlockSize.Resize(numbl);
	int64_t ibl,iblcount,graphsize = 0;
	for(ibl=0, iblcount=0; ibl<numbl; ibl++)
	{
		if(colors[ibl]==color) 
		{
			this->fBlockSize[iblcount++] = blockgraphindex[ibl+1]-blockgraphindex[ibl];
			graphsize += this->fBlockSize[iblcount-1];
		}
	}
	fBlock.Resize(graphsize);
	fBlockIndex.Resize(iblcount+1);
	fBlockIndex[0] = 0;
	for(ibl=0, iblcount=0; ibl<numbl; ibl++)
	{
		if(colors[ibl]==color) 
		{
			int64_t first = blockgraphindex[ibl];
			int64_t last = blockgraphindex[ibl+1];
			int64_t firstcp = fBlockIndex[iblcount];
			this->fBlockIndex[iblcount+1] = firstcp+this->fBlockSize[iblcount];
			//      int lastcp = fBlockIndex[iblcount+1];
			int64_t ieq,ieqcp;
			for(ieq=first,ieqcp=firstcp; ieq<last; ieq++,ieqcp++)
			{
				fBlock[ieqcp] = blockgraph[ieq];
			}
			iblcount++;
		}
	}
	this->fBlockSize.Resize(iblcount);
	this->Initialize(this->fBlockSize);
	this->fRow = rows;
	this->fCol = rows;
}


template<class TVar>
TPZSparseBlockDiagonal<TVar>::~TPZSparseBlockDiagonal()
{
}

template<class TVar>
const TVar TPZSparseBlockDiagonal<TVar>::Get(const int64_t row, const int64_t col) const
{
	int64_t rblock,rblockindex,cblock,cblockindex;
	FindBlockIndex(row,rblock,rblockindex);
	if(rblock == -1) return (TVar)0;
	FindBlockIndex(col,cblock,cblockindex);
	if(cblock != rblock) return (TVar)0;
	int64_t pos = rblockindex + cblockindex*this->fBlockSize[rblock];
	return this->fStorage[this->fBlockPos[rblock]+pos];
}

template<class TVar>
const TVar TPZSparseBlockDiagonal<TVar>::GetVal(const int64_t row, const int64_t col) const
{
	int64_t rblock,rblockindex,cblock,cblockindex;
	FindBlockIndex(row,rblock,rblockindex);
	if(rblock == -1) return (TVar) 0;
	FindBlockIndex(col,cblock,cblockindex);
	if(cblock != rblock) return (TVar) 0;
	int64_t pos = rblockindex + cblockindex*this->fBlockSize[rblock];
	return this->fStorage[this->fBlockPos[rblock]+pos];
}

template <class TVar>
int TPZSparseBlockDiagonal<TVar>::Put(const int64_t row, const int64_t col, const TVar& value)
{
	int64_t rblock,rblockindex,cblock,cblockindex;
	FindBlockIndex(row,rblock,rblockindex);
	if(rblock == -1) return -1;
	FindBlockIndex(col,cblock,cblockindex);
	if(cblock != rblock) return -1;
	int64_t pos = rblockindex + cblockindex*this->fBlockSize[rblock];
	this->fStorage[this->fBlockPos[rblock]+pos] = value;
	return 0;
}

template<class TVar>
int TPZSparseBlockDiagonal<TVar>::PutVal(const int64_t row, const int64_t col, const TVar& value)
{
	int64_t rblock,rblockindex,cblock,cblockindex;
	FindBlockIndex(row,rblock,rblockindex);
	if(rblock == -1) return -1;
	FindBlockIndex(col,cblock,cblockindex);
	if(cblock != rblock) return -1;
	int64_t pos = rblockindex + cblockindex*this->fBlockSize[rblock];
	this->fStorage[this->fBlockPos[rblock]+pos] = value;
	return 0;
}
template<class TVar>
TVar& TPZSparseBlockDiagonal<TVar>::operator ( )(const int64_t row, const int64_t col)
{
	int64_t rblock,rblockindex,cblock,cblockindex;
	FindBlockIndex(row,rblock,rblockindex);
	if(rblock == -1) return this->gZero;
	FindBlockIndex(col,cblock,cblockindex);
	if(cblock != rblock) return this->gZero;
	int64_t pos = rblockindex + cblockindex*this->fBlockSize[rblock];
	return this->fStorage[this->fBlockPos[rblock]+pos];
}

template<class TVar>
int TPZSparseBlockDiagonal<TVar>::Substitution(TPZFMatrix<TVar>* B) const
{
	TPZFNMatrix<1000,TVar > BG(fBlock.NElements(),B->Cols());
	Gather(*B,BG);
	int result = TPZBlockDiagonal<TVar>::Substitution(&BG);
	B->Zero();
	Scatter(BG,*B);
	return result;
	
}

template<class TVar>
TVar& TPZSparseBlockDiagonal<TVar>::s(const int64_t row, const int64_t col)
{
	int64_t rblock,rblockindex,cblock,cblockindex;
	FindBlockIndex(row,rblock,rblockindex);
	if(rblock == -1) return this->gZero;
	FindBlockIndex(col,cblock,cblockindex);
	if(cblock != rblock) return this->gZero;
	int64_t pos = rblockindex + cblockindex*this->fBlockSize[rblock];
	return this->fStorage[this->fBlockPos[rblock]+pos];
}

template<class TVar>
void TPZSparseBlockDiagonal<TVar>::Print(const char* message, std::ostream& out, const MatrixOutputFormat format) const
{
	TPZBlockDiagonal<TVar>::Print(message, out, format);
	if(format == EFormatted)
	{
		out << "Equations for each block " << endl;
		int64_t nbl = fBlockIndex.NElements()-1;
		int64_t ibl;
		for(ibl = 0; ibl<nbl ; ibl++)
		{
			int64_t first = fBlockIndex[ibl];
			int64_t last = fBlockIndex[ibl+1];
			out << "Block " << ibl << " : ";
			int64_t i;
			for(i=first; i<last; i++) out << fBlock[i] << " ";
			out << endl;
		}
	}
}

template<class TVar>
void TPZSparseBlockDiagonal<TVar>::AddBlock(int64_t i, TPZFMatrix<TVar>& block)
{
    TPZBlockDiagonal<TVar>::AddBlock(i, block);
}

template<class TVar>
void TPZSparseBlockDiagonal<TVar>::BuildFromMatrix(TPZMatrix<TVar>& matrix)
{
#ifdef PZ_LOG
	LOGPZ_DEBUG(logger, "TPZSparseBlockDiagonal::BuildFromMatrix");
#endif
	TPZManVector<int64_t> indices;
	TPZFNMatrix<10000,TVar> submat(0,0);
	int64_t ibl,nbl = fBlockIndex.NElements()-1;
	for(ibl=0; ibl<nbl; ibl++)
	{
		int64_t nel = this->fBlockSize[ibl];
		indices.Resize(nel);
		submat.Resize(nel,nel);
		int64_t iel,first = fBlockIndex[ibl];
		for(iel=0; iel<nel; iel++) indices[iel] = fBlock[first+iel];
		matrix.GetSub(indices,submat);
		this->SetBlock(ibl,submat);
	}
}

template<class TVar>
void TPZSparseBlockDiagonal<TVar>::GetBlock(int64_t i, TPZFMatrix<TVar>& block)
{
    TPZBlockDiagonal<TVar>::GetBlock(i, block);
}

template<class TVar>
void TPZSparseBlockDiagonal<TVar>::MultAdd(const TPZFMatrix<TVar>& x, const TPZFMatrix<TVar>& y, TPZFMatrix<TVar>& z, const TVar alpha, const TVar beta, const int opt) const
{
#ifdef PZ_LOG
	LOGPZ_DEBUG(logger, "TPZSparseBlockDiagonal::MultAdd");
#endif
	TPZFNMatrix<1000,TVar> xsc(0,0),ysc(0,0,0.),zsc(0,0);
	xsc.Resize(this->fBlock.NElements(),x.Cols());
	z.Zero();
	if(abs(beta) != 0.) ysc.Resize(fBlock.NElements(),y.Cols());
	zsc.Resize(fBlock.NElements(),z.Cols());
	Gather(x,xsc);
	if(abs(beta) != 0.) Scatter(y,ysc);
	TPZBlockDiagonal<TVar>::MultAdd(xsc, ysc, zsc, alpha, beta, opt);
	Scatter(zsc,z);
}

/*!
 \fn TPZSparseBlockDiagonal::FindBlockIndex(int glob, int &block, int &blockind)
 */
template<class TVar>
void TPZSparseBlockDiagonal<TVar>::FindBlockIndex(int64_t glob, int64_t &block, int64_t &blockind) const
{
    int64_t numbl = fBlockIndex.NElements()-2;
    int64_t ieq,ibl;
    for(ibl = 0; ibl<numbl; ibl++)
    {
		for(ieq = fBlockIndex[ibl];ieq<fBlockIndex[ibl+1];ieq++)
		{
			if(fBlock[ieq] == glob)
			{
				block = ibl;
				blockind = ieq-fBlockIndex[ibl];
				return;
			}
		}
    }
    block = -1;
    blockind = -1;
}

/*!
 \fn TPZSparseBlockDiagonal::Scatter(TPZFMatrix<>&in, TPZFMatrix<>&out) const
 */
template<class TVar>
void TPZSparseBlockDiagonal<TVar>::Scatter(const TPZFMatrix<TVar> &in, TPZFMatrix<TVar> &out) const
{
    int64_t neq = fBlock.NElements();
    int64_t nc = in.Cols();
    int64_t ieq,ic;
    for(ic=0; ic<nc; ic++)
    {
		for(ieq=0; ieq<neq; ieq++) out(fBlock[ieq],ic) += in.GetVal(ieq,ic);
    }
}

/*!
 \fn TPZSparseBlockDiagonal::Gather(TPZFMatrix<>&in, TPZFMatrix<>&out) const
 */
template<class TVar>
void TPZSparseBlockDiagonal<TVar>::Gather(const TPZFMatrix<TVar> &in, TPZFMatrix<TVar> &out) const
{
    int64_t neq = fBlock.NElements();
    int64_t nc = in.Cols();
    int64_t ieq,ic;
    for(ic=0; ic<nc; ic++)
    {
		for(ieq=0; ieq<neq; ieq++) out(ieq,ic) = in.GetVal(fBlock[ieq],ic);
    }
}

/**
 * Updates the values of the matrix based on the values of the matrix
 */
template<class TVar>
void TPZSparseBlockDiagonal<TVar>::UpdateFrom(TPZAutoPointer<TPZMatrix<TVar> > mat)
{
#ifdef PZ_LOG
	LOGPZ_DEBUG(logger, "TPZSparseBlockDiagonal::UpdateFrom");
#endif
	if(!mat) 
	{
		cout << __PRETTY_FUNCTION__ << " called with zero argument\n";
		return;
	}
	this->fDecomposed = ENoDecompose;
	int64_t nblock = this->fBlockSize.NElements();
	int64_t b,bsize,pos;
	TPZManVector<int64_t,1000> indices;
	for(b=0; b<nblock; b++) {
		bsize = this->fBlockSize[b];
		indices.Resize(bsize);
		int64_t r;
		pos = this->fBlockPos[b];
		for(r=0; r<bsize; r++) indices[r] = fBlock[fBlockIndex[b]+r]; 
		TPZFMatrix<TVar> block(bsize,bsize,&this->fStorage[pos],bsize*bsize);
		mat->GetSub(indices,block);
	}
	
}

template <class TVar>
int TPZSparseBlockDiagonal<TVar>::ClassId() const{
    return Hash("TPZSparseBlockDiagonal") ^ TPZBlockDiagonal<TVar>::ClassId() << 1;
}
template class TPZSparseBlockDiagonal<float>;
template class TPZSparseBlockDiagonal<double>;
template class TPZSparseBlockDiagonal<long double>;

template class TPZSparseBlockDiagonal<std::complex<float> >;
template class TPZSparseBlockDiagonal<std::complex<double> >;
template class TPZSparseBlockDiagonal<std::complex<long double> >;
