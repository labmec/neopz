/**
 * @file
 * @brief Contains the implementation of the TPZSparseBlockDiagonal methods.
 */
//
// C++ Implementation: tpzsparseblockdiagonal
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2004
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "tpzsparseblockdiagonal.h"
#include "pzfmatrix.h"

#include "pzlog.h"
#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.StrMatrix"));
#endif

using namespace std;

TPZSparseBlockDiagonal::TPZSparseBlockDiagonal()
{
}

TPZSparseBlockDiagonal::TPZSparseBlockDiagonal(TPZVec<int> &blockgraph, TPZVec<int> &blockgraphindex,int rows) : fBlock(blockgraph), fBlockIndex(blockgraphindex)
{
	int numbl = blockgraphindex.NElements()-1;
	fBlockSize.Resize(numbl);
	int ibl;
	for(ibl=0; ibl<numbl; ibl++)
	{
		fBlockSize[ibl] = blockgraphindex[ibl+1]-blockgraphindex[ibl];
	}
	Initialize(fBlockSize);
	fRow = rows;
	fCol = rows;
}

TPZSparseBlockDiagonal::TPZSparseBlockDiagonal(TPZVec<int> &blockgraph, TPZVec<int> &blockgraphindex,int rows, int color, TPZVec<int> &colors)
{
	LOGPZ_DEBUG(logger, "Constructor of TPZSparseBlockDiagonal");
	int numbl = blockgraphindex.NElements()-1;
	fBlockSize.Resize(numbl);
	int ibl,iblcount,graphsize = 0;
	for(ibl=0, iblcount=0; ibl<numbl; ibl++)
	{
		if(colors[ibl]==color) 
		{
			fBlockSize[iblcount++] = blockgraphindex[ibl+1]-blockgraphindex[ibl];
			graphsize += fBlockSize[iblcount-1];
		}
	}
	fBlock.Resize(graphsize);
	fBlockIndex.Resize(iblcount+1);
	fBlockIndex[0] = 0;
	for(ibl=0, iblcount=0; ibl<numbl; ibl++)
	{
		if(colors[ibl]==color) 
		{
			int first = blockgraphindex[ibl];
			int last = blockgraphindex[ibl+1];
			int firstcp = fBlockIndex[iblcount];
			fBlockIndex[iblcount+1] = firstcp+fBlockSize[iblcount];
			//      int lastcp = fBlockIndex[iblcount+1];
			int ieq,ieqcp;
			for(ieq=first,ieqcp=firstcp; ieq<last; ieq++,ieqcp++)
			{
				fBlock[ieqcp] = blockgraph[ieq];
			}
			iblcount++;
		}
	}
	fBlockSize.Resize(iblcount);
	Initialize(fBlockSize);
	fRow = rows;
	fCol = rows;
}



TPZSparseBlockDiagonal::~TPZSparseBlockDiagonal()
{
}

//static REAL gZero = 0.;

const REAL& TPZSparseBlockDiagonal::Get(const int row, const int col) const
{
	int rblock,rblockindex,cblock,cblockindex;
	FindBlockIndex(row,rblock,rblockindex);
	if(rblock == -1) return gZero;
	FindBlockIndex(col,cblock,cblockindex);
	if(cblock != rblock) return gZero;
	int pos = rblockindex + cblockindex*fBlockSize[rblock];
	return fStorage[fBlockPos[rblock]+pos];
}

const REAL& TPZSparseBlockDiagonal::GetVal(const int row, const int col) const
{
	int rblock,rblockindex,cblock,cblockindex;
	FindBlockIndex(row,rblock,rblockindex);
	if(rblock == -1) return gZero;
	FindBlockIndex(col,cblock,cblockindex);
	if(cblock != rblock) return gZero;
	int pos = rblockindex + cblockindex*fBlockSize[rblock];
	return fStorage[fBlockPos[rblock]+pos];
}

int TPZSparseBlockDiagonal::Put(const int row, const int col, const REAL& value)
{
	int rblock,rblockindex,cblock,cblockindex;
	FindBlockIndex(row,rblock,rblockindex);
	if(rblock == -1) return -1;
	FindBlockIndex(col,cblock,cblockindex);
	if(cblock != rblock) return -1;
	int pos = rblockindex + cblockindex*fBlockSize[rblock];
	fStorage[fBlockPos[rblock]+pos] = value;
	return 0;
}

int TPZSparseBlockDiagonal::PutVal(const int row, const int col, const REAL& value)
{
	int rblock,rblockindex,cblock,cblockindex;
	FindBlockIndex(row,rblock,rblockindex);
	if(rblock == -1) return -1;
	FindBlockIndex(col,cblock,cblockindex);
	if(cblock != rblock) return -1;
	int pos = rblockindex + cblockindex*fBlockSize[rblock];
	fStorage[fBlockPos[rblock]+pos] = value;
	return 0;
}

REAL& TPZSparseBlockDiagonal::operator ( )(const int row, const int col)
{
	int rblock,rblockindex,cblock,cblockindex;
	FindBlockIndex(row,rblock,rblockindex);
	if(rblock == -1) return gZero;
	FindBlockIndex(col,cblock,cblockindex);
	if(cblock != rblock) return gZero;
	int pos = rblockindex + cblockindex*fBlockSize[rblock];
	return fStorage[fBlockPos[rblock]+pos];
}

int TPZSparseBlockDiagonal::Substitution(TPZFMatrix* B) const
{
	TPZFNMatrix<1000> BG(fBlock.NElements(),B->Cols());
	Gather(*B,BG,1);
	int result = TPZBlockDiagonal::Substitution(&BG);
	B->Zero();
	Scatter(BG,*B,1);
	return result;
	
}

REAL& TPZSparseBlockDiagonal::s(const int row, const int col)
{
	int rblock,rblockindex,cblock,cblockindex;
	FindBlockIndex(row,rblock,rblockindex);
	if(rblock == -1) return gZero;
	FindBlockIndex(col,cblock,cblockindex);
	if(cblock != rblock) return gZero;
	int pos = rblockindex + cblockindex*fBlockSize[rblock];
	return fStorage[fBlockPos[rblock]+pos];
}

void TPZSparseBlockDiagonal::Print(const char* message, std::ostream& out, const MatrixOutputFormat format) const
{
	TPZBlockDiagonal::Print(message, out, format);
	if(format == EFormatted)
	{
		out << "Equations for each block " << endl;
		int nbl = fBlockIndex.NElements()-1;
		int ibl;
		for(ibl = 0; ibl<nbl ; ibl++)
		{
			int first = fBlockIndex[ibl];
			int last = fBlockIndex[ibl+1];
			out << "Block " << ibl << " : ";
			int i;
			for(i=first; i<last; i++) out << fBlock[i] << " ";
			out << endl;
		}
	}
}

void TPZSparseBlockDiagonal::AddBlock(int i, TPZFMatrix& block)
{
    TPZBlockDiagonal::AddBlock(i, block);
}

void TPZSparseBlockDiagonal::BuildFromMatrix(TPZMatrix& matrix)
{
	LOGPZ_DEBUG(logger, "TPZSparseBlockDiagonal::BuildFromMatrix");
	TPZManVector<int> indices;
	TPZFNMatrix<10000> submat(0,0);
	int ibl,nbl = fBlockIndex.NElements()-1;
	for(ibl=0; ibl<nbl; ibl++)
	{
		int nel = fBlockSize[ibl];
		indices.Resize(nel);
		submat.Resize(nel,nel);
		int iel,first = fBlockIndex[ibl];
		for(iel=0; iel<nel; iel++) indices[iel] = fBlock[first+iel];
		matrix.GetSub(indices,submat);
		SetBlock(ibl,submat);
	}
}

void TPZSparseBlockDiagonal::GetBlock(int i, TPZFMatrix& block)
{
    TPZBlockDiagonal::GetBlock(i, block);
}

void TPZSparseBlockDiagonal::MultAdd(const TPZFMatrix& x, const TPZFMatrix& y, TPZFMatrix& z, const REAL alpha, const REAL beta, const int opt, const int stride) const
{
	LOGPZ_DEBUG(logger, "TPZSparseBlockDiagonal::MultAdd");
	TPZFNMatrix<1000000> xsc(0,0),ysc(0,0,0.),zsc(0,0);
	xsc.Resize(fBlock.NElements(),x.Cols());
	z.Zero();
	if(beta != 0.) ysc.Resize(fBlock.NElements(),y.Cols());
	zsc.Resize(fBlock.NElements(),z.Cols());
	//  int fRowKeep = fRow;
	//  int fColKeep = fCol;
	//  fRow = fBlock.NElements();
	//  fCol = fBlock.NElements();
	Gather(x,xsc,stride);
	if(beta != 0.) Scatter(y,ysc,stride);
	TPZBlockDiagonal::MultAdd(xsc, ysc, zsc, alpha, beta, opt,1);
	Scatter(zsc,z,stride);
	//  fRow = fRowKeep;
	//  fCol = fColKeep;
}



/*!
 \fn TPZSparseBlockDiagonal::FindBlockIndex(int glob, int &block, int &blockind)
 */
void TPZSparseBlockDiagonal::FindBlockIndex(int glob, int &block, int &blockind) const
{
    int numbl = fBlockIndex.NElements()-2;
    int ieq,ibl;
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
 \fn TPZSparseBlockDiagonal::Scatter(TPZFMatrix &in, TPZFMatrix &out) const
 */
void TPZSparseBlockDiagonal::Scatter(const TPZFMatrix &in, TPZFMatrix &out, int stride) const
{
    int neq = fBlock.NElements();
    int nc = in.Cols();
    int ieq,ic;
    for(ic=0; ic<nc; ic++)
    {
		for(ieq=0; ieq<neq; ieq++) out(fBlock[ieq]*stride,ic) += in.GetVal(ieq,ic);
    }
}


/*!
 \fn TPZSparseBlockDiagonal::Gather(TPZFMatrix &in, TPZFMatrix &out, int stride) const
 */
void TPZSparseBlockDiagonal::Gather(const TPZFMatrix &in, TPZFMatrix &out, int stride) const
{
    int neq = fBlock.NElements();
    int nc = in.Cols();
    int ieq,ic;
    for(ic=0; ic<nc; ic++)
    {
		for(ieq=0; ieq<neq; ieq++) out(ieq,ic) = in.GetVal(fBlock[ieq]*stride,ic);
    }
}

/**
 * Updates the values of the matrix based on the values of the matrix
 */
void TPZSparseBlockDiagonal::UpdateFrom(TPZAutoPointer<TPZMatrix> mat)
{
	LOGPZ_DEBUG(logger, "TPZSparseBlockDiagonal::UpdateFrom");
	if(!mat) 
	{
		cout << __PRETTY_FUNCTION__ << " called with zero argument\n";
		return;
	}
	this->fDecomposed = ENoDecompose;
	int nblock = fBlockSize.NElements();
	int b,bsize,pos;
	TPZManVector<int,1000> indices;
	for(b=0; b<nblock; b++) {
		bsize = fBlockSize[b];
		indices.Resize(bsize);
		int r;
		pos = fBlockPos[b];
		for(r=0; r<bsize; r++) indices[r] = fBlock[fBlockIndex[b]+r]; 
		TPZFMatrix block(bsize,bsize,&fStorage[pos],bsize*bsize);
		mat->GetSub(indices,block);
	}
	
}

