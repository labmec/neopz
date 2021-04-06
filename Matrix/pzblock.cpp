/**
 * @file
 * @brief Contains the implementation of the TPZBlock methods. Purpose: Let to see selected block of a matrix.
 */

#include <stdlib.h>
#include <stdio.h>
#include "pzerror.h"
#include "pzblock.h"
#include "TPZStream.h"

#include <sstream>
#include "pzlog.h"
#ifdef PZ_LOG
static TPZLogger logger("pz.matrix.tpzblock");
#endif

using namespace std;
template<class TVar>
REAL TPZBlock<TVar>::gZero = 0;


/*************************** Public ***************************/

/******************/
/*** Construtor ***/
template<class TVar>
TPZBlock<TVar>::TPZBlock( TPZMatrix<TVar> *const pMatrix,const int nBlocks,const int dim )
: TPZRegisterClassId(&TPZBlock::ClassId) 
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
TPZBlock<TVar>::TPZBlock(const TPZBlock<TVar> &bl) : TPZRegisterClassId(&TPZBlock::ClassId), 
fBlock(bl.fBlock) {
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
	const auto maxBlocks = fBlock.NElements();
	if ( (bRow >= maxBlocks) || (bCol >= maxBlocks) ){
		PZError << __PRETTY_FUNCTION__ <<'\n';
		PZError << "Get <block index out of range>" <<std::endl;
		DebugStop();
	}
	
	const auto rowDim = fBlock[bRow].dim;
	const auto colDim = fBlock[bCol].dim;
	
	
	if ( !rowDim || !colDim ){
		PZError << __PRETTY_FUNCTION__ << '\n';
		PZError << "Get <inexistent block>";
		PZError << std::endl;
		DebugStop();
	}
	
	if ( (r >= rowDim) || (c >= colDim) ){
		PZError << __PRETTY_FUNCTION__<<'\n';
		PZError << "Get <element is out of the block>" <<std::endl;
		DebugStop();
	}
	return this->GetVal(bRow,bCol,r,c);
}

/***********/
/*** Put ***/
template<class TVar>
int
TPZBlock<TVar>::Put(const int bRow,const int bCol,const int r,const int c,
					const TVar& value )
{
	const auto maxBlocks = fBlock.NElements();
	if ( (bRow >= maxBlocks) || (bCol >= maxBlocks) ){
		PZError << __PRETTY_FUNCTION__ <<'\n';
		PZError << "Get <block index out of range>" <<std::endl;
		DebugStop();
	}
	
	const auto rowDim = fBlock[bRow].dim;
	const auto colDim = fBlock[bCol].dim;
	
	
	if ( !rowDim || !colDim ){
		PZError << __PRETTY_FUNCTION__ << '\n';
		PZError << "Put <inexistent block>";
		PZError << std::endl;
		DebugStop();
	}
	
	if ( (r >= rowDim) || (c >= colDim) ){
		PZError << __PRETTY_FUNCTION__<<'\n';
		PZError << "Put <element is out of the block>" <<std::endl;
		DebugStop();
	}
	return this->PutVal(bRow, bCol, r, c, value);
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
		PZError<<__PRETTY_FUNCTION__<< "Get <block index out of range>" <<std::endl; DebugStop();
	}
	
	int rowDim = fBlock[bRow].dim;
	
	
	if ( !rowDim  )
	{
		PZError<<__PRETTY_FUNCTION__<< "Get <inexistent block>" <<std::endl; DebugStop();
	}
	
	if ( (row >= rowDim) || (col >= fpMatrix->Cols()) )
	{
		PZError<<__PRETTY_FUNCTION__<< "Get <elemente is out of the block>" <<std::endl; DebugStop();
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
		{PZError<<__PRETTY_FUNCTION__<< "Put <block index out of range>" <<std::endl; DebugStop();}
	
	int rowDim = fBlock[bRow].dim;
	
	if ( !rowDim )
		{PZError<<__PRETTY_FUNCTION__<< "Put <inexistent block>" <<std::endl; DebugStop();}
	
	if ( (row >= rowDim) || (col >= fpMatrix->Cols()) )
		{PZError<<__PRETTY_FUNCTION__<< "Put <elemente is out of the block>" <<std::endl; DebugStop();}
	
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


template<class TVar>
void TPZBlock<TVar>::PrintStructure(std::ostream &out)
{
    int MaxBlocks = fBlock.NElements();
    for ( int bRow = 0; bRow < MaxBlocks; bRow++ )
    {
        out << "row block " << bRow << " pos " << fBlock[bRow].pos << " dim " << fBlock[bRow].dim << "\n";
    }
}


/*************************** Private ***************************/

template<class TVar>
void TPZBlock<TVar>::TNode::Read(TPZStream &buf, void *context) { //ok
    buf.Read(&pos,1);
    buf.Read(&dim,1);
}

template<class TVar>
void TPZBlock<TVar>::TNode::Write(TPZStream &buf, int withclassid) const { //ok
    buf.Write(&pos,1);
    buf.Write(&dim,1);
}

#ifndef BORLAND
template class TPZRestoreClass< TPZBlock<float> >;
template class TPZRestoreClass< TPZBlock<double> >;
template class TPZRestoreClass< TPZBlock<long double> >;

template class TPZRestoreClass< TPZBlock<std::complex<float> > >;
template class TPZRestoreClass< TPZBlock<std::complex<double> > >;
template class TPZRestoreClass< TPZBlock<std::complex<long double> > >;
#endif


/** Saves the element data to a stream */
template<class TVar>
void TPZBlock<TVar>::Write(TPZStream &buf, int withclassid) const { //ok
	buf.Write(fBlock);
    TPZPersistenceManager::WritePointer(fpMatrix, &buf);
	
}

/** Reads the element data from a stream */
template<class TVar>
void TPZBlock<TVar>::Read(TPZStream &buf, void *context) { //ok
	buf.Read<TNode>(fBlock,context);
    fpMatrix = dynamic_cast<TPZMatrix<TVar> *>(TPZPersistenceManager::GetInstance(&buf));
}

template class TPZBlock<float>;
template class TPZBlock<double>;
template class TPZBlock<long double>;

template class TPZBlock<std::complex<float> >;
template class TPZBlock<std::complex<double> >;
template class TPZBlock<std::complex<long double> >;
