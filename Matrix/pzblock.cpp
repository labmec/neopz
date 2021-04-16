/**
 * @file
 * @brief Contains the implementation of the TPZBlock methods. Purpose: Let to see selected block of a matrix.
 */

#include <stdlib.h>
#include <stdio.h>
#include "pzerror.h"
#include "pzblock.h"
#include "TPZStream.h"
#include "pzfmatrix.h"
#include <sstream>
#include "pzlog.h"
#ifdef PZ_LOG
static TPZLogger logger("pz.matrix.tpzblock");
#endif

using namespace std;
REAL TPZBlock::gZero = 0;


/*************************** Public ***************************/

/******************/
/*** Construtor ***/

TPZBlock::TPZBlock( TPZBaseMatrix *const pMatrix,const int nBlocks,const int dim )
: TPZRegisterClassId(&TPZBlock::ClassId) 
{
	const int maxBlocks = [&]() -> int {
        if(pMatrix) return  nBlocks ? nBlocks : pMatrix->Rows();
        else return nBlocks;
    }();
	
	//  fBlock = 0;
	if(maxBlocks) fBlock.Resize(maxBlocks);
	fpMatrix = pMatrix;
	

	const int mat_size = [&]() -> int {
        if(pMatrix) {return pMatrix->Rows();}
        else {return 1;}
    }();
    
	const int dim2 = [&]() ->int {
        // The row dimension of the matrix determines the size of the block object
        if(pMatrix && dim*nBlocks!=mat_size){return mat_size/maxBlocks;}
        else{return dim;}
    }();

	int pos = 0;
	for ( int i = 0; i < fBlock.NElements(); i++, pos += dim2 )
    {
		fBlock[i].pos = pos;
		fBlock[i].dim = dim2;
    }
	if(maxBlocks && dim2) fBlock[maxBlocks-1].dim = dim2 + mat_size%dim2;
	else if(maxBlocks) fBlock[maxBlocks-1].dim = mat_size-dim2;
}


TPZBlock::TPZBlock(const TPZBlock &bl) : TPZRegisterClassId(&TPZBlock::ClassId), 
fBlock(bl.fBlock) {
	fpMatrix = bl.fpMatrix;
}

/*******************/
/*** operator =  ***/

TPZBlock &TPZBlock::operator=(const TPZBlock & bl) {
	if(this == &bl) return *this;
	fBlock = bl.fBlock;
	fpMatrix = bl.fpMatrix;
	return *this;
}

/******************/
/***  Destrutor ***/

TPZBlock::~TPZBlock() {
}


/******************/
/*** Set Blocks ***/

int
TPZBlock::SetNBlocks(const int num_of_blocks )
{
	//modified Philippe 24/7/97
	// a small optimization
	int MaxBlocks = fBlock.NAlloc();
	if(num_of_blocks >= MaxBlocks) fBlock.Expand((int) (num_of_blocks*1.2));
	TNode copy;
	fBlock.Resize(num_of_blocks,copy);
	if(num_of_blocks == MaxBlocks) return 1;
	return ( 1 );
}


int TPZBlock::Set(const int b,const int dim,const int pos ) {
	if ( b >= fBlock.NElements() ) {
		cout << "TPZBlock::Set called with parameter out of range\n";
		return( 0 );
	}
	
	fBlock[b].dim = dim;
	if ( pos >= 0 )
		fBlock[b].pos = pos;
	else
		fBlock[b].pos = (b ? fBlock[b-1].pos + fBlock[b-1].dim : 0);
	
	return( 1 );
}


/**************/
/*** SetAll ***/

int
TPZBlock::SetAll( TPZVec<int> & dimensions )
{
	int total_dim=0;
	int i,nel = dimensions.NElements() ;
	for(i=0;i<nel;i++) total_dim += dimensions[i];
	if ( total_dim != fpMatrix->Rows() ||
		 total_dim > fpMatrix->Rows() ){
		PZError<<__PRETTY_FUNCTION__<<"SetAll <new block dimensions not compatible whit matrix dimension>"<<std::endl;
		DebugStop();
	}
	
	int pos=0;
	
	//SetNumBlocks( nel );
	SetNBlocks( nel );
	
	for(i=0; i<nel;i++ )
    {
		fBlock[i].pos=pos;
		pos=pos+(fBlock[i].dim=dimensions[i]);
    }
	
	return ( 1 );
}

/*****************/
/*** Resequence **/

int TPZBlock::Resequence(const int start) {
	int MaxBlocks = fBlock.NElements();
	if (start>=MaxBlocks) return 0;
	for (int i= start+1; i < MaxBlocks; i++)
		fBlock[i].pos=fBlock[i-1].pos+fBlock[i-1].dim;
	return ( 1 );
}

/**************/
/*** Remove ***/

int
TPZBlock::Remove(const int index )
{
	int MaxBlocks = fBlock.NElements();
	if ( index >= MaxBlocks )
		return( 0 );
	
	fBlock[index].dim = 0;    // and the corresponding elements into the fpMatrix???
	return( 1 );
}

/***************/
/*** Verify ****/

int
TPZBlock::Verify() const
{
	int MaxBlocks = fBlock.NElements();
	for ( int i = 0; i < MaxBlocks-1; i++ )
		if (fBlock[i].pos + fBlock[i].dim != fBlock[i+1].pos) return ( 0 );
	
	if (fBlock[MaxBlocks-1].pos + fBlock[MaxBlocks-1].dim != fpMatrix->Rows())
		return ( 0 );
	return ( 1 );
}



/// Return the index in the blocked matrix

int64_t TPZBlock::Index(const int64_t bRow, const int r) const
{
    auto MaxBlocks = fBlock.NElements();
    int64_t row(r);
    if(bRow <0 || bRow >= MaxBlocks || row < 0 || row >= fBlock[bRow].dim) {
        PZError << __PRETTY_FUNCTION__ <<" indexes out of range\n";
        DebugStop();
    }
    row += fBlock[bRow].pos;
    return row;
}


void TPZBlock::PrintStructure(std::ostream &out)
{
    int MaxBlocks = fBlock.NElements();
    for ( int bRow = 0; bRow < MaxBlocks; bRow++ )
    {
        out << "row block " << bRow << " pos " << fBlock[bRow].pos << " dim " << fBlock[bRow].dim << "\n";
    }
}


template<class TVar>
int
TPZBlock::PutBlock(const int bRow,const int bCol,const TPZFMatrix<TVar> & block )
{
    auto tmp = this->Matrix<TVar>();
	return( tmp->PutSub( fBlock[bRow].pos,
							 fBlock[bCol].pos, block ) );
}

/*****************/
/*** Get Block ***/
template<class TVar>
int
TPZBlock::GetBlock(const int bRow,const int bCol, TPZFMatrix<TVar> &block ) const
{
	const auto row = fBlock[bRow].pos;
	const auto col = fBlock[bCol].pos;
	const auto rowDim = fBlock[bRow].dim;
	const auto colDim = fBlock[bCol].dim;
    auto tmp = this->Matrix<TVar>();
	if ( rowDim && colDim )
		return( tmp->GetSub( row, col, rowDim, colDim, block ) );
	else
		return( 0 );
}

/*************************** Private ***************************/


void TPZBlock::TNode::Read(TPZStream &buf, void *context) { //ok
    buf.Read(&pos,1);
    buf.Read(&dim,1);
}


void TPZBlock::TNode::Write(TPZStream &buf, int withclassid) const { //ok
    buf.Write(&pos,1);
    buf.Write(&dim,1);
}

#ifndef BORLAND
template class TPZRestoreClass<TPZBlock>;
#endif


/** Saves the element data to a stream */

void TPZBlock::Write(TPZStream &buf, int withclassid) const { //ok
	buf.Write(fBlock);
    TPZPersistenceManager::WritePointer(fpMatrix, &buf);
	
}

/** Reads the element data from a stream */

void TPZBlock::Read(TPZStream &buf, void *context) { //ok
	buf.Read<TNode>(fBlock,context);
    fpMatrix = dynamic_cast<TPZBaseMatrix *>(TPZPersistenceManager::GetInstance(&buf));
}

int TPZBlock::ClassId() const {
    return Hash("TPZBlock");
}

template TPZFMatrix<float> * TPZBlock::Matrix();
template TPZFMatrix<double> * TPZBlock::Matrix();
template TPZFMatrix<long double> * TPZBlock::Matrix();

template TPZFMatrix<std::complex<float> > * TPZBlock::Matrix();
template TPZFMatrix<std::complex<double> > * TPZBlock::Matrix();
template TPZFMatrix<std::complex<long double> > * TPZBlock::Matrix();

template const TPZFMatrix<float> * TPZBlock::Matrix() const;
template const TPZFMatrix<double> * TPZBlock::Matrix() const;
template const TPZFMatrix<long double> * TPZBlock::Matrix() const;

template const TPZFMatrix<std::complex<float> > * TPZBlock::Matrix() const;
template const TPZFMatrix<std::complex<double> > * TPZBlock::Matrix() const;
template const TPZFMatrix<std::complex<long double> > * TPZBlock::Matrix() const;



template
int TPZBlock::PutBlock(const int, const int,
                       const TPZFMatrix<float> &);
template
int TPZBlock::PutBlock(const int, const int,
                       const TPZFMatrix<double> &);
template
int TPZBlock::PutBlock(const int, const int,
                       const TPZFMatrix<long double> &);
template
int TPZBlock::PutBlock(const int, const int,
                       const TPZFMatrix<std::complex<float>> &);
template
int TPZBlock::PutBlock(const int, const int,
                       const TPZFMatrix<std::complex<double>> &);
template
int TPZBlock::PutBlock(const int, const int,
                       const TPZFMatrix<std::complex<long double>> &);


template
int TPZBlock::GetBlock(const int, const int,
                       TPZFMatrix<float> &) const;
template
int TPZBlock::GetBlock(const int, const int,
                       TPZFMatrix<double> &) const;
template
int TPZBlock::GetBlock(const int, const int,
                       TPZFMatrix<long double> &) const;
template
int TPZBlock::GetBlock(const int, const int,
                       TPZFMatrix<std::complex<float>> &) const;
template
int TPZBlock::GetBlock(const int, const int,
                       TPZFMatrix<std::complex<double>> &) const;
template
int TPZBlock::GetBlock(const int, const int,
                       TPZFMatrix<std::complex<long double>> &) const;
