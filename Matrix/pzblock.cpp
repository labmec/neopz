
//
// Author: MISAEL LUIS SANTANA MANDUJANO.
//
// File:   tblock.cc
//
// Class:  TPZBlock
//
// Obs.:   Permite a visualizacao de matrizes atraves de blocos.
//
// Versao: 12 / 1994.

#include <stdlib.h>
#include <stdio.h>
#include "pzerror.h"
#include "pzblock.h"

REAL TPZBlock::gZero = 0;//Cedric


/*************************** Public ***************************/

/******************/
/*** Construtor ***/

TPZBlock::TPZBlock( TPZMatrix *const pMatrix,const int nBlocks,const int dim )
{
  if(pMatrix) fMaxBlocks = ( nBlocks ? nBlocks : pMatrix->Rows() );
  else fMaxBlocks = nBlocks;

  fBlock = 0;
  if(fMaxBlocks) fBlock = new( TNode[fMaxBlocks] );
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
      dim2 = mat_size/fMaxBlocks;
  }
  // fim de sugestao

  int pos = 0;
  for ( int i = 0; i < fMaxBlocks; i++, pos += dim2 )
    {
      fBlock[i].pos = pos;
      fBlock[i].dim = dim2;
    }
  if(fMaxBlocks && dim2) fBlock[fMaxBlocks-1].dim = dim2 + mat_size%dim2;
  else if(fMaxBlocks) fBlock[fMaxBlocks-1].dim = mat_size-dim2;
}

TPZBlock::TPZBlock(const TPZBlock &bl) {
  fMaxBlocks = bl.fMaxBlocks;
  fBlock = new( TNode[fMaxBlocks] );
  int ibl;
  for(ibl=0; ibl<fMaxBlocks; ibl++) {
    fBlock[ibl].pos = bl.fBlock[ibl].pos;
    fBlock[ibl].dim = bl.fBlock[ibl].dim;
  }
  fpMatrix = bl.fpMatrix;
}

/*******************/
/*** operator =  ***/
TPZBlock &TPZBlock::operator=(const TPZBlock & bl) {
  if(this == &bl) return *this;
  if(fMaxBlocks != bl.fMaxBlocks) {
    delete [] fBlock;
    fBlock = new( TNode[bl.fMaxBlocks] );
  }
  fMaxBlocks = bl.fMaxBlocks;
  int ibl;
  for(ibl=0; ibl<fMaxBlocks; ibl++) {
    fBlock[ibl] = bl.fBlock[ibl];
  }
  fpMatrix = bl.fpMatrix;
  return *this;
}

/******************/
/***  Destrutor ***/
TPZBlock::~TPZBlock() {
  delete[] fBlock;
}

/******************/
/*** SetMatrix ****/
void
TPZBlock::SetMatrix(TPZMatrix *const other){
  fpMatrix = other;
}


/******************/
/*** Set Blocks ***/
int
TPZBlock::SetNBlocks(const int num_of_blocks )
{
  //modified Philippe 24/7/97
  // a small optimization
  if(num_of_blocks == fMaxBlocks) return 1;
  if ( num_of_blocks <= fMaxBlocks )
    {
      int i;
      for (i = num_of_blocks; i < fMaxBlocks; i++ ){
	fBlock[i].pos=0;fBlock[i].dim=0;
      }
      fMaxBlocks = num_of_blocks;
      return( 1 );
    }

  TNode *newBlocks = new( TNode[num_of_blocks] );
  int i;
  int lastpos = Dim();
  for ( i = 0; i < fMaxBlocks; i++ )
    newBlocks[i] = fBlock[i];
  for ( ; i < num_of_blocks; i++ ) {
    newBlocks[i].dim = 0;
    newBlocks[i].pos = lastpos;
  }

  delete []fBlock;
  fBlock     = newBlocks;
  fMaxBlocks = num_of_blocks;
  return ( 1 );
}

int TPZBlock::Set(const int b,const int dim,const int pos ) {
  if ( b >= fMaxBlocks ) {
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
       total_dim > fpMatrix->Rows() )
    Error("SetAll <new block dimensions not compatible whit matrix dimension>");

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
  if (start>=fMaxBlocks) return 0;
  for (int i= start+1; i < fMaxBlocks; i++)
    fBlock[i].pos=fBlock[i-1].pos+fBlock[i-1].dim;
  return ( 1 );
}

/**************/
/*** Remove ***/
int
TPZBlock::Remove(const int index )
{
  if ( index >= fMaxBlocks )
    return( 0 );

  fBlock[index].dim = 0;    // and the corresponding elements into the fpMatrix???
  return( 1 );
}



/***************/
/*** Verify ****/
int
TPZBlock::Verify() const
{
  for ( int i = 0; i < fMaxBlocks-1; i++ )
    if (fBlock[i].pos + fBlock[i].dim != fBlock[i+1].pos) return ( 0 );

  if (fBlock[fMaxBlocks-1].pos + fBlock[fMaxBlocks-1].dim != fpMatrix->Rows())
    return ( 0 );
  return ( 1 );
}


/***********/
/*** Get ***/
const REAL &
TPZBlock::Get(const int bRow,const int bCol,const int r,const int c ) const
{
  int row(r),col(c);
  if ( (bRow >= fMaxBlocks) || (bCol >= fMaxBlocks) )
    Error( "Get <block index out of range>" );

  int rowDim = fBlock[bRow].dim;
  int colDim = fBlock[bCol].dim;


  if ( !rowDim || !colDim )
    Error( "Get <inexistent block>" );

  if ( (row >= rowDim) || (col >= colDim) )
    Error( "Get <elemente is out of the block>" );

  row += fBlock[bRow].pos;
  col += fBlock[bCol].pos;
  return( fpMatrix->Get( row, col ) );
}



/***********/
/*** Put ***/
int
TPZBlock::Put(const int bRow,const int bCol,const int r,const int c,
	      const REAL& value )
{
  int row(r),col(c);
  if ( (bRow >= fMaxBlocks) || (bCol >= fMaxBlocks) )
    Error( "Put <block index out of range>" );

  int rowDim = fBlock[bRow].dim;
  int colDim = fBlock[bCol].dim;

  if ( !rowDim || !colDim )
    Error( "Put <inexistent block>" );

  if ( (row >= rowDim) || (col >= colDim) )
    Error( "Put <elemente is out of the block>" );

  row += fBlock[bRow].pos;
  col += fBlock[bCol].pos;
  return( fpMatrix->Put( row, col, value ) );
}



/**************/
/*** GetVal ***/
const REAL &
TPZBlock::GetVal(const int bRow,const int bCol,const int r,const int c ) const
{
  int row(r),col(c);
  if(bRow <0 || bRow >= fMaxBlocks || bCol <0 || bCol >= fMaxBlocks || row < 0 || row >= fBlock[bRow].dim) {
    cout << "TPZBlock::GetVal indexes out of range\n";
    exit(-1);
  }
  row += fBlock[bRow].pos;
  col += fBlock[bCol].pos;
  return( fpMatrix->Get( row, col ) );
}


REAL &
TPZBlock::operator()(const int bRow,const int bCol,const int r,const int c )
{
  int row(r),col(c);
  if(bRow <0 || bRow >= fMaxBlocks || bCol <0 || bCol >= fMaxBlocks || row < 0 || row >= fBlock[bRow].dim) {
    cout << "TPZBlock::operator() indexes out of range\n";
    exit(-1);
  }
  row += fBlock[bRow].pos;
  col += fBlock[bCol].pos;
  return( (*fpMatrix)( row, col ) );
}



/**************/
/*** PutVal ***/
int
TPZBlock::PutVal(const int bRow,const int bCol,const int r,const int c,
		 const REAL& value )
{
  int row(r),col(c);
  row += fBlock[bRow].pos;
  col += fBlock[bCol].pos;
  return( fpMatrix->Put( row, col, value ) );
}



/*****************/
/*** Put Block ***/
int
TPZBlock::PutBlock(const int bRow,const int bCol,const TPZFMatrix & block )
{
  return( fpMatrix->PutSub( fBlock[bRow].pos,
			    fBlock[bCol].pos, block ) );
}



/*****************/
/*** Get Block ***/
int
TPZBlock::GetBlock(const int bRow,const int bCol, TPZFMatrix *const block ) const
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
int
TPZBlock::AddBlock(const int bRow,const int bCol,const TPZFMatrix& block )
{
  return( fpMatrix->AddSub( fBlock[bRow].pos,
			    fBlock[bCol].pos, block ) );
}
/********************/
/**** InsertBLock****/
int
TPZBlock::InsertBlock(const int block_row,const int block_col,
		      const int target_row,const int target_col, TPZMatrix &target) const
{
  int rowDim = fBlock[block_row].dim;
  int colDim = fBlock[block_col].dim;
  int row = target_row;
  int col = target_col;

  if ( ((target_row + rowDim) > target.Rows()) ||
       ((target_col + colDim) > target.Cols()) ) {
    Error( "GetSub <the sub-matrix is too big>" ) ;
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
int
TPZBlock::PrintBlock(const int bRow,const int bCol,const char *title,
		     TPZostream &out ) const
{
  out << title << ":";

  for ( int r = 0; r < fBlock[bRow].dim; r++ )
    {
      out << "\n  ";
      for ( int c = 0; c < fBlock[bCol].dim; c++ )
	out << GetVal( bRow, bCol, r, c ) << "  ";
    }
  out << "\n";
  return( 1 );
}



/*************/
/*** Print ***/
void
TPZBlock::Print(const char *title, TPZostream &out,TPZMatrix *mat) {
  TPZMatrix *sol=fpMatrix;
  if (mat) SetMatrix(mat);
  char block_title[32];

  out << title << ":";
  for ( int bRow = 0; bRow < fMaxBlocks; bRow++ )
    for ( int bCol = 0; bCol < fMaxBlocks; bCol++ )
      {
	out << "\n";
	sprintf( block_title, "Block (%d,%d) of %dX%d:", bRow, bCol,
		 fBlock[bRow].dim,fBlock[bCol].dim );
	PrintBlock(bRow,bCol,block_title,out);
      }
  out << "\n";
  SetMatrix( sol);
}



/*************************** Private ***************************/

/*************/
/*** Error ***/
void
TPZBlock::Error(const char *msg ) const
{
  cout << "TPZBlock::" << msg << ".\n";
  // pzerror.Show();
  exit( 1 );
}

#ifdef OOPARLIB

int TPZBlock::Unpack (TReceiveStorage *buf ){
  TSaveable::Unpack(buf);
  fBlock= new TNode [fMaxBlocks];
  buf->UpkInt((int *)fBlock,2*fMaxBlocks);
  TSaveable *sav = buf->Restore();
  if(!sav->DerivedFrom("TPZMatrix")) exit(-1);
  fpMatrix = (TPZMatrix *) sav;
  return 1;
}


TSaveable *TPZBlock::Restore(TReceiveStorage *buf) {
  TPZBlock *m = new TPZBlock(0,0,0);
  m->Unpack(buf);
  return m;
}

int TPZBlock::Pack( TSendStorage *buf ) const {
  TSaveable::Pack(buf);
  buf->PkInt( (int *) fBlock,2*fMaxBlocks);
  fpMatrix->Pack(buf);
  return 1;
}


int TPZBlock::DerivedFrom(const long Classid) const {
  if(Classid == GetClassID()) return 1;
  return TSaveable::DerivedFrom(Classid);
}

int TPZBlock::DerivedFrom(const char *classname) const {

  if(!strcmp(ClassName(),classname)) return 1;
  return TSaveable::DerivedFrom(classname);
}

#endif


