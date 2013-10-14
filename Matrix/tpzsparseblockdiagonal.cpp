
/**
 * @file
 * @brief Contains the implementation of the TPZSparseBlockDiagonal methods.
 */

#include "tpzsparseblockdiagonal.h"
#include "pzfmatrix.h"

#include "pzlog.h"
#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.StrMatrix"));
#endif

using namespace std;

template<class TVar>
TPZSparseBlockDiagonal<TVar>::TPZSparseBlockDiagonal()
{
}
template<class TVar>
TPZSparseBlockDiagonal<TVar>::TPZSparseBlockDiagonal(TPZVec<long> &blockgraph, TPZVec<long> &blockgraphindex,long rows) : fBlock(blockgraph), fBlockIndex(blockgraphindex)
{
	long numbl = blockgraphindex.NElements()-1;
	this->fBlockSize.Resize(numbl);
	long ibl;
	for(ibl=0; ibl<numbl; ibl++)
	{
		this->fBlockSize[ibl] = blockgraphindex[ibl+1]-blockgraphindex[ibl];
	}
	this->Initialize(this->fBlockSize);
	this->fRow = rows;
	this->fCol = rows;
}

template<class TVar>
TPZSparseBlockDiagonal<TVar>::TPZSparseBlockDiagonal(TPZVec<long> &blockgraph, TPZVec<long> &blockgraphindex,long rows, int color, TPZVec<int> &colors)
{
	LOGPZ_DEBUG(logger, "Constructor of TPZSparseBlockDiagonal");
	long numbl = blockgraphindex.NElements()-1;
	this->fBlockSize.Resize(numbl);
	long ibl,iblcount,graphsize = 0;
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
			long first = blockgraphindex[ibl];
			long last = blockgraphindex[ibl+1];
			long firstcp = fBlockIndex[iblcount];
			this->fBlockIndex[iblcount+1] = firstcp+this->fBlockSize[iblcount];
			//      int lastcp = fBlockIndex[iblcount+1];
			long ieq,ieqcp;
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
const TVar& TPZSparseBlockDiagonal<TVar>::Get(const long row, const long col) const
{
	long rblock,rblockindex,cblock,cblockindex;
	FindBlockIndex(row,rblock,rblockindex);
	if(rblock == -1) return this->gZero;
	FindBlockIndex(col,cblock,cblockindex);
	if(cblock != rblock) return this->gZero;
	long pos = rblockindex + cblockindex*this->fBlockSize[rblock];
	return this->fStorage[this->fBlockPos[rblock]+pos];
}

template<class TVar>
const TVar& TPZSparseBlockDiagonal<TVar>::GetVal(const long row, const long col) const
{
	long rblock,rblockindex,cblock,cblockindex;
	FindBlockIndex(row,rblock,rblockindex);
	if(rblock == -1) return this->gZero;
	FindBlockIndex(col,cblock,cblockindex);
	if(cblock != rblock) return this->gZero;
	long pos = rblockindex + cblockindex*this->fBlockSize[rblock];
	return this->fStorage[this->fBlockPos[rblock]+pos];
}

template <class TVar>
int TPZSparseBlockDiagonal<TVar>::Put(const long row, const long col, const TVar& value)
{
	long rblock,rblockindex,cblock,cblockindex;
	FindBlockIndex(row,rblock,rblockindex);
	if(rblock == -1) return -1;
	FindBlockIndex(col,cblock,cblockindex);
	if(cblock != rblock) return -1;
	long pos = rblockindex + cblockindex*this->fBlockSize[rblock];
	this->fStorage[this->fBlockPos[rblock]+pos] = value;
	return 0;
}

template<class TVar>
int TPZSparseBlockDiagonal<TVar>::PutVal(const long row, const long col, const TVar& value)
{
	long rblock,rblockindex,cblock,cblockindex;
	FindBlockIndex(row,rblock,rblockindex);
	if(rblock == -1) return -1;
	FindBlockIndex(col,cblock,cblockindex);
	if(cblock != rblock) return -1;
	long pos = rblockindex + cblockindex*this->fBlockSize[rblock];
	this->fStorage[this->fBlockPos[rblock]+pos] = value;
	return 0;
}
template<class TVar>
TVar& TPZSparseBlockDiagonal<TVar>::operator ( )(const long row, const long col)
{
	long rblock,rblockindex,cblock,cblockindex;
	FindBlockIndex(row,rblock,rblockindex);
	if(rblock == -1) return this->gZero;
	FindBlockIndex(col,cblock,cblockindex);
	if(cblock != rblock) return this->gZero;
	long pos = rblockindex + cblockindex*this->fBlockSize[rblock];
	return this->fStorage[this->fBlockPos[rblock]+pos];
}

template<class TVar>
int TPZSparseBlockDiagonal<TVar>::Substitution(TPZFMatrix<TVar>* B) const
{
	TPZFNMatrix<1000,TVar > BG(fBlock.NElements(),B->Cols());
	Gather(*B,BG,1);
	int result = TPZBlockDiagonal<TVar>::Substitution(&BG);
	B->Zero();
	Scatter(BG,*B,1);
	return result;
	
}

template<class TVar>
TVar& TPZSparseBlockDiagonal<TVar>::s(const long row, const long col)
{
	long rblock,rblockindex,cblock,cblockindex;
	FindBlockIndex(row,rblock,rblockindex);
	if(rblock == -1) return this->gZero;
	FindBlockIndex(col,cblock,cblockindex);
	if(cblock != rblock) return this->gZero;
	long pos = rblockindex + cblockindex*this->fBlockSize[rblock];
	return this->fStorage[this->fBlockPos[rblock]+pos];
}

template<class TVar>
void TPZSparseBlockDiagonal<TVar>::Print(const char* message, std::ostream& out, const MatrixOutputFormat format) const
{
	TPZBlockDiagonal<TVar>::Print(message, out, format);
	if(format == EFormatted)
	{
		out << "Equations for each block " << endl;
		long nbl = fBlockIndex.NElements()-1;
		long ibl;
		for(ibl = 0; ibl<nbl ; ibl++)
		{
			long first = fBlockIndex[ibl];
			long last = fBlockIndex[ibl+1];
			out << "Block " << ibl << " : ";
			long i;
			for(i=first; i<last; i++) out << fBlock[i] << " ";
			out << endl;
		}
	}
}

template<class TVar>
void TPZSparseBlockDiagonal<TVar>::AddBlock(long i, TPZFMatrix<TVar>& block)
{
    TPZBlockDiagonal<TVar>::AddBlock(i, block);
}

template<class TVar>
void TPZSparseBlockDiagonal<TVar>::BuildFromMatrix(TPZMatrix<TVar>& matrix)
{
	LOGPZ_DEBUG(logger, "TPZSparseBlockDiagonal::BuildFromMatrix");
	TPZManVector<long> indices;
	TPZFNMatrix<10000,TVar> submat(0,0);
	long ibl,nbl = fBlockIndex.NElements()-1;
	for(ibl=0; ibl<nbl; ibl++)
	{
		long nel = this->fBlockSize[ibl];
		indices.Resize(nel);
		submat.Resize(nel,nel);
		long iel,first = fBlockIndex[ibl];
		for(iel=0; iel<nel; iel++) indices[iel] = fBlock[first+iel];
		matrix.GetSub(indices,submat);
		this->SetBlock(ibl,submat);
	}
}

template<class TVar>
void TPZSparseBlockDiagonal<TVar>::GetBlock(long i, TPZFMatrix<TVar>& block)
{
    TPZBlockDiagonal<TVar>::GetBlock(i, block);
}

template<class TVar>
void TPZSparseBlockDiagonal<TVar>::MultAdd(const TPZFMatrix<TVar>& x, const TPZFMatrix<TVar>& y, TPZFMatrix<TVar>& z, const TVar alpha, const TVar beta, const int opt, const int stride) const
{
	LOGPZ_DEBUG(logger, "TPZSparseBlockDiagonal::MultAdd");
	TPZFNMatrix<1000,TVar> xsc(0,0),ysc(0,0,0.),zsc(0,0);
	xsc.Resize(this->fBlock.NElements(),x.Cols());
	z.Zero();
	if(abs(beta) != 0.) ysc.Resize(fBlock.NElements(),y.Cols());
	zsc.Resize(fBlock.NElements(),z.Cols());
	Gather(x,xsc,stride);
	if(abs(beta) != 0.) Scatter(y,ysc,stride);
	TPZBlockDiagonal<TVar>::MultAdd(xsc, ysc, zsc, alpha, beta, opt,1);
	Scatter(zsc,z,stride);
}

/*!
 \fn TPZSparseBlockDiagonal::FindBlockIndex(int glob, int &block, int &blockind)
 */
template<class TVar>
void TPZSparseBlockDiagonal<TVar>::FindBlockIndex(long glob, long &block, long &blockind) const
{
    long numbl = fBlockIndex.NElements()-2;
    long ieq,ibl;
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
void TPZSparseBlockDiagonal<TVar>::Scatter(const TPZFMatrix<TVar> &in, TPZFMatrix<TVar> &out, int stride) const
{
    long neq = fBlock.NElements();
    long nc = in.Cols();
    long ieq,ic;
    for(ic=0; ic<nc; ic++)
    {
		for(ieq=0; ieq<neq; ieq++) out(fBlock[ieq]*stride,ic) += in.GetVal(ieq,ic);
    }
}

/*!
 \fn TPZSparseBlockDiagonal::Gather(TPZFMatrix<>&in, TPZFMatrix<>&out, int stride) const
 */
template<class TVar>
void TPZSparseBlockDiagonal<TVar>::Gather(const TPZFMatrix<TVar> &in, TPZFMatrix<TVar> &out, int stride) const
{
    long neq = fBlock.NElements();
    long nc = in.Cols();
    long ieq,ic;
    for(ic=0; ic<nc; ic++)
    {
		for(ieq=0; ieq<neq; ieq++) out(ieq,ic) = in.GetVal(fBlock[ieq]*stride,ic);
    }
}

/**
 * Updates the values of the matrix based on the values of the matrix
 */
template<class TVar>
void TPZSparseBlockDiagonal<TVar>::UpdateFrom(TPZAutoPointer<TPZMatrix<TVar> > mat)
{
	LOGPZ_DEBUG(logger, "TPZSparseBlockDiagonal::UpdateFrom");
	if(!mat) 
	{
		cout << __PRETTY_FUNCTION__ << " called with zero argument\n";
		return;
	}
	this->fDecomposed = ENoDecompose;
	long nblock = this->fBlockSize.NElements();
	long b,bsize,pos;
	TPZManVector<long,1000> indices;
	for(b=0; b<nblock; b++) {
		bsize = this->fBlockSize[b];
		indices.Resize(bsize);
		long r;
		pos = this->fBlockPos[b];
		for(r=0; r<bsize; r++) indices[r] = fBlock[fBlockIndex[b]+r]; 
		TPZFMatrix<TVar> block(bsize,bsize,&this->fStorage[pos],bsize*bsize);
		mat->GetSub(indices,block);
	}
	
}

template class TPZSparseBlockDiagonal<float>;
template class TPZSparseBlockDiagonal<double>;
template class TPZSparseBlockDiagonal<long double>;

template class TPZSparseBlockDiagonal<std::complex<float> >;
template class TPZSparseBlockDiagonal<std::complex<double> >;
template class TPZSparseBlockDiagonal<std::complex<long double> >;
