/**
 * @file
 * @brief Contains the implementation of the TPZBlock methods. Purpose: Let to see selected block of a matrix.
 */

#include <stdlib.h>
#include <stdio.h>
#include "pzerror.h"
#include "pzblock.h"
#include "pzstream.h"

#include <sstream>
#include "pzlog.h"
#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.matrix.tpzblock"));
#endif

using namespace std;
template<class TVar>
REAL TPZBlock<TVar>::gZero = 0;


/*************************** Public ***************************/

/******************/
/*** Construtor ***/
template<class TVar>
TPZBlock<TVar>::TPZBlock( TPZMatrix<TVar> *const pMatrix,const int nBlocks,const int dim )
{
	int MaxBlocks = 0;
	if(pMatrix) MaxBlocks = ( nBlocks ? nBlocks : pMatrix->Rows() );
	else MaxBlocks = nBlocks;
	
	//  fBlock = 0;
	if(MaxBlocks) fBlock.Resize(MaxBlocks);
	fpMatrix = pMatrix;
	
	// sugestao de implementacao:
	// o construtor adaptaria a dimensao das matrizes da diagonal ao
	// numero de blocos: por exemplo, 4 blocos em uma matriz M 8x8
	// ( TPZBlock(&M,4) ) teriam dimensao 2.
	//  PROBLEMA: cria uma variavel interna, dim2
	int dim2 = dim;
	int mat_size = 1;
	if(pMatrix) {
		// The row dimension of the matrix determines the size of the block object
		mat_size = pMatrix->Rows();
		if ( (dim*nBlocks!=mat_size) )
			dim2 = mat_size/MaxBlocks;
	}
	// fim de sugestao
	
	int pos = 0;
	for ( int i = 0; i < fBlock.NElements(); i++, pos += dim2 )
    {
		fBlock[i].pos = pos;
		fBlock[i].dim = dim2;
    }
	if(MaxBlocks && dim2) fBlock[MaxBlocks-1].dim = dim2 + mat_size%dim2;
	else if(MaxBlocks) fBlock[MaxBlocks-1].dim = mat_size-dim2;
}

template<class TVar>
TPZBlock<TVar>::TPZBlock(const TPZBlock<TVar> &bl) : fBlock(bl.fBlock) {
	fpMatrix = bl.fpMatrix;
}

/*******************/
/*** operator =  ***/
template<class TVar>
TPZBlock<TVar> &TPZBlock<TVar>::operator=(const TPZBlock<TVar> & bl) {
	if(this == &bl) return *this;
	fBlock = bl.fBlock;
	fpMatrix = bl.fpMatrix;
	return *this;
}

/******************/
/***  Destrutor ***/
template<class TVar>
TPZBlock<TVar>::~TPZBlock() {
}

/******************/
/*** SetMatrix ****/
template<class TVar>
void
TPZBlock<TVar>::SetMatrix(TPZMatrix<TVar> *const other){
	fpMatrix = other;
}


/******************/
/*** Set Blocks ***/
template<class TVar>
int
TPZBlock<TVar>::SetNBlocks(const int num_of_blocks )
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

template<class TVar>
int TPZBlock<TVar>::Set(const int b,const int dim,const int pos ) {
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
template<class TVar>
int
TPZBlock<TVar>::SetAll( TPZVec<int> & dimensions )
{
	int total_dim=0;
	int i,nel = dimensions.NElements() ;
	for(i=0;i<nel;i++) total_dim += dimensions[i];
	if ( total_dim != fpMatrix->Rows() ||
		total_dim > fpMatrix->Rows() )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"SetAll <new block dimensions not compatible whit matrix dimension>");
	
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
template<class TVar>
int TPZBlock<TVar>::Resequence(const int start) {
	int MaxBlocks = fBlock.NElements();
	if (start>=MaxBlocks) return 0;
	for (int i= start+1; i < MaxBlocks; i++)
		fBlock[i].pos=fBlock[i-1].pos+fBlock[i-1].dim;
	return ( 1 );
}

/**************/
/*** Remove ***/
template<class TVar>
int
TPZBlock<TVar>::Remove(const int index )
{
	int MaxBlocks = fBlock.NElements();
	if ( index >= MaxBlocks )
		return( 0 );
	
	fBlock[index].dim = 0;    // and the corresponding elements into the fpMatrix???
	return( 1 );
}

/***************/
/*** Verify ****/
template<class TVar>
int
TPZBlock<TVar>::Verify() const
{
	int MaxBlocks = fBlock.NElements();
	for ( int i = 0; i < MaxBlocks-1; i++ )
		if (fBlock[i].pos + fBlock[i].dim != fBlock[i+1].pos) return ( 0 );
	
	if (fBlock[MaxBlocks-1].pos + fBlock[MaxBlocks-1].dim != fpMatrix->Rows())
		return ( 0 );
	return ( 1 );
}

/***********/
/*** Get ***/
template<class TVar>
const TVar &
TPZBlock<TVar>::Get(const int bRow,const int bCol,const int r,const int c ) const
{
	int row(r),col(c);
	int MaxBlocks = fBlock.NElements();
	if ( (bRow >= MaxBlocks) || (bCol >= MaxBlocks) )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Get <block index out of range>" );
	
	int rowDim = fBlock[bRow].dim;
	int colDim = fBlock[bCol].dim;
	
	
	if ( !rowDim || !colDim )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Get <inexistent block>" );
	
	if ( (row >= rowDim) || (col >= colDim) )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Get <elemente is out of the block>" );
	
	row += fBlock[bRow].pos;
	col += fBlock[bCol].pos;
	return( fpMatrix->Get( row, col ) );
}

/***********/
/*** Put ***/
template<class TVar>
int
TPZBlock<TVar>::Put(const int bRow,const int bCol,const int r,const int c,
					const TVar& value )
{
	int MaxBlocks = fBlock.NElements();
	int row(r),col(c);
	if ( (bRow >= MaxBlocks) || (bCol >= MaxBlocks) )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Put <block index out of range>" );
	
	int rowDim = fBlock[bRow].dim;
	int colDim = fBlock[bCol].dim;
	
	if ( !rowDim || !colDim )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Put <inexistent block>" );
	
	if ( (row >= rowDim) || (col >= colDim) )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Put <elemente is out of the block>" );
	
	row += fBlock[bRow].pos;
	col += fBlock[bCol].pos;
	return( fpMatrix->Put( row, col, value ) );
}

/***********/
/*** Get ***/
template<class TVar>
const TVar &
TPZBlock<TVar>::Get(const int bRow,const int r,const int c ) const
{
	int row(r),col(c);
	int MaxBlocks = fBlock.NElements();
	if ( (bRow >= MaxBlocks)  )
	{
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Get <block index out of range>" );
	}
	
	int rowDim = fBlock[bRow].dim;
	
	
	if ( !rowDim  )
	{
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Get <inexistent block>" );
	}
	
	if ( (row >= rowDim) || (col >= fpMatrix->Cols()) )
	{
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Get <elemente is out of the block>" );
	}
	
	row += fBlock[bRow].pos;
	return( fpMatrix->Get( row, col ) );
}

/***********/
/*** Put ***/
template<class TVar>
int
TPZBlock<TVar>::Put(const int bRow,const int r,const int c,
					const TVar& value )
{
	int MaxBlocks = fBlock.NElements();
	int row(r),col(c);
	if ( (bRow >= MaxBlocks)  )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Put <block index out of range>" );
	
	int rowDim = fBlock[bRow].dim;
	
	if ( !rowDim )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Put <inexistent block>" );
	
	if ( (row >= rowDim) || (col >= fpMatrix->Cols()) )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Put <elemente is out of the block>" );
	
	row += fBlock[bRow].pos;
	return( fpMatrix->Put( row, col, value ) );
}

/**************/
/*** GetVal ***/
template<class TVar>
const TVar &
TPZBlock<TVar>::GetVal(const int bRow,const int bCol,const int r,const int c ) const
{
	int MaxBlocks = fBlock.NElements();
	int row(r),col(c);
	if(bRow <0 || bRow >= MaxBlocks || bCol <0 || bCol >= MaxBlocks || row < 0 || row >= fBlock[bRow].dim) {
		cout << "TPZBlock::GetVal indexes out of range\n";
		DebugStop();
	}
	row += fBlock[bRow].pos;
	col += fBlock[bCol].pos;
	return( fpMatrix->Get( row, col ) );
}

template<class TVar>
TVar &
TPZBlock<TVar>::operator()(const int bRow,const int bCol,const int r,const int c ) const
{
	int MaxBlocks = fBlock.NElements();
	int row(r),col(c);
	if(bRow <0 || bRow >= MaxBlocks || bCol <0 || bCol >= MaxBlocks || row < 0 || row >= fBlock[bRow].dim) {
		cout << "TPZBlock::operator() indexes out of range\n";
		DebugStop();
	}
	row += fBlock[bRow].pos;
	col += fBlock[bCol].pos;
	return( (*fpMatrix)( row, col ) );
}

/**************/
/*** PutVal ***/
template<class TVar>
int
TPZBlock<TVar>::PutVal(const int bRow,const int bCol,const int r,const int c,
					   const TVar& value )
{
	int row(r),col(c);
	row += fBlock[bRow].pos;
	col += fBlock[bCol].pos;
	return( fpMatrix->Put( row, col, value ) );
}

/*****************/
/*** Put Block ***/
template<class TVar>
int
TPZBlock<TVar>::PutBlock(const int bRow,const int bCol,const TPZFMatrix<TVar> & block )
{
	return( fpMatrix->PutSub( fBlock[bRow].pos,
							 fBlock[bCol].pos, block ) );
}

/*****************/
/*** Get Block ***/
template<class TVar>
int
TPZBlock<TVar>::GetBlock(const int bRow,const int bCol, TPZFMatrix<TVar> *const block ) const
{
	int row = fBlock[bRow].pos;
	int col = fBlock[bCol].pos;
	int rowDim = fBlock[bRow].dim;
	int colDim = fBlock[bCol].dim;
	if ( rowDim && colDim )
		return( fpMatrix->GetSub( row, col, rowDim, colDim, *block ) );
	else
		return( 0 );
}

/*****************/
/*** Add Block ***/
template<class TVar>
int
TPZBlock<TVar>::AddBlock(const int bRow,const int bCol,const TPZFMatrix<TVar>& block )
{
	return( fpMatrix->AddSub( fBlock[bRow].pos,
							 fBlock[bCol].pos, block ) );
}

/********************/
/**** InsertBLock****/
template<class TVar>
int
TPZBlock<TVar>::InsertBlock(const int block_row,const int block_col,
							const int target_row,const int target_col, TPZMatrix<TVar> &target) const
{
	int rowDim = fBlock[block_row].dim;
	int colDim = fBlock[block_col].dim;
	int row = target_row;
	int col = target_col;
	
	if ( ((target_row + rowDim) > target.Rows()) ||
		((target_col + colDim) > target.Cols()) ) {
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "GetSub <the sub-matrix is too big>" ) ;
		return ( 0 );
	}
	
	for ( int r = 0; r < rowDim; r++,row++){
		int pcol=col;
		for ( int c = 0; c < colDim; c++,pcol++ )
			target.PutVal( row, pcol, GetVal(block_row , block_col , r, c ) );
	}
	
	return( 1 );
}

/*******************/
/*** Print Block ***/
template<class TVar>
int
TPZBlock<TVar>::PrintBlock(const int bRow,const int bCol,const char *title,
						   TPZostream &out ) const
{
	out << title << ":";
	
	for ( int r = 0; r < fBlock[bRow].dim; r++ )
    {
		out << "\n  ";
		for ( int c = 0; c < 1/*fBlock[bCol].dim*/; c++ )
			out << GetVal( bRow, bCol, r, c ) << "  ";
    }
	out << "\n";
	return( 1 );
}

/*************/
/*** Print ***/
template<class TVar>
void
TPZBlock<TVar>::Print(const char *title, TPZostream &out,TPZMatrix<TVar> *mat) {
	TPZMatrix<TVar> *sol=fpMatrix;
	if (mat) SetMatrix(mat);
	char block_title[32];
	
	int MaxBlocks = fBlock.NElements();
	out << title << ":\n";
	for ( int bRow = 0; bRow < MaxBlocks; bRow++ )
	{
		out << "block " << bRow << " pos " << fBlock[bRow].pos << " dim " << fBlock[bRow].dim << "\n";  
		for ( int bCol = 0; bCol < 1/*MaxBlocks*/; bCol++ )
		{
			out << "\n";
			sprintf( block_title, "Block (%d,%d) of %dX%d:", bRow, bCol,
					fBlock[bRow].dim,fBlock[bCol].dim );
			PrintBlock(bRow,bCol,block_title,out);
		}
	}
	out << "\n";
	SetMatrix( sol);
}

/*************/
/*** Print ***/
template<class TVar>
void
TPZBlock<TVar>::PrintSolution(const char *title, TPZostream &out) {
	TPZMatrix<TVar> *sol=fpMatrix;
	
	char block_title[32];
	
	int MaxBlocks = fBlock.NElements();
	out << title << ":";
	for ( int bRow = 0; bRow < MaxBlocks; bRow++ )
	{
		out << "\n";
		sprintf( block_title, "Block (%d,%d) of %dX%d:", bRow, 0,
				fBlock[bRow].dim,fBlock[0].dim );
		PrintBlock(bRow,0,block_title,out);
	}
	out << "\n";
	SetMatrix( sol);
}

/*************************** Private ***************************/

template<class TVar>
void TPZBlock<TVar>::TNode::Read(TPZStream &buf, void *context)
{
	buf.Read(&dim,1);
	buf.Read(&pos,1);
}

template<class TVar>
void TPZBlock<TVar>::TNode::Write(TPZStream &buf, void *context)
{
	buf.Write(&dim,1);
	buf.Write(&pos,1);
}

/** Returns the unique identifier for reading/writing objects to streams */
template<class TVar>
int TPZBlock<TVar>::ClassId() const
{
	return TPZBLOCKID;
}


template class TPZRestoreClass< TPZBlock<REAL> , TPZBLOCKID>;

/** Saves the element data to a stream */
template<class TVar>
void TPZBlock<TVar>::Write(TPZStream &buf, int withclassid)
{
	TPZSaveable::Write(buf,withclassid);
	TPZSaveable::WriteObjects<TNode>(buf,fBlock);
	
}

/** Reads the element data from a stream */
template<class TVar>
void TPZBlock<TVar>::Read(TPZStream &buf, void *context)
{
	fpMatrix = (TPZMatrix<TVar> *) context;
	TPZSaveable::Read(buf,context);
	ReadObjects<TNode>(buf,fBlock,context);
}

template class TPZBlock<REAL>;
