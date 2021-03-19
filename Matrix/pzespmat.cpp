/**
 * @file
 * @brief Contains the implementation of the TPZSpMatrix methods.
 */

#include <math.h>
#include "pzespmat.h"
#include "pzfmatrix.h"

#include "pzerror.h"
#include "stdlib.h"
#include "pzlink.h"

#include <sstream>
#include "pzlog.h"
#ifdef PZ_LOG
static PZLogger logger("pz.matrix.tpzspmatrix");
#endif

using namespace std;

/*******************/
/*** TPZSpMatrix ***/

/**************************** PUBLIC ****************************/
/**************************/
/*** Construtor (int) ***/
template<class TVar>
TPZSpMatrix<TVar>::TPZSpMatrix(const int64_t rows,const int64_t cols )
: TPZRegisterClassId(&TPZSpMatrix::ClassId),
TPZMatrix<TVar>( rows, cols )
#ifdef WORKPOOL
, fWp()
#endif
{
	fElem = new TPZLink<TPZNode>[ rows ] ;
	if ( fElem == NULL )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "TPZSpMatrix( dim ) <Error creating Matrix>" );
#ifdef WORKPOOL
	for(int64_t i=0; i<rows; i++) fElem[i].SetWorkPool(&fWp);
#endif
}


/*****************/
/*** Destrutor ***/

template<class TVar>
TPZSpMatrix<TVar>::~TPZSpMatrix ()
{
	Clear();
}


/***********/
/*** Put ***/

template<class TVar>
int
TPZSpMatrix<TVar>::Put(const int64_t row,const int64_t col,const TVar& value )
{
	if ( (row >= this->Rows()) || (col >= this->Cols()) || row <0 || col<0)
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Put <indices out of band matrix range>" );
	
	return( PutVal( row, col, value ) );
}



/***********/
/*** Get ***/

template<class TVar>
const TVar &
TPZSpMatrix<TVar>::Get(const int64_t row,const int64_t col ) const
{
	if ( (row >= this->Rows()) || (col >= this->Cols()) || row<0 || col<0)
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Get <indices out of band matrix range>" );
	
	return( GetVal( row, col ) );
}



/**************/
/*** PutVal ***/
//
//  Escreve um elemento na matriz, sem verificar fronteiras.
//  O valor da linha (row) deve ser maior ou igual ao da coluna
// (col).
//
template<class TVar>
int
TPZSpMatrix<TVar>::PutVal(const int64_t row,const int64_t col,const TVar & value )
{
	TPZLink<TPZNode> *pRow = &fElem[row];
	TPZNode        node;
	
	// Caso a lista esteja vazia, garante que (node.col != col).
	node.col = -1;
	
	// Procura pela posicao que esta ou deveria estar o elemento.
	pRow->Head();
	while ( pRow->Get( &node ) && (node.col < col) )
		pRow->Next();
	
	node.elem = value;
	
	// Se encontrou a posicao do elemento...
	if ( node.col == col )
    {
		// Se o elemento e' ZERO, remove-o.
		if ( IsZero( value ) )
			pRow->Remove();
		
		// Se o elemento nao e' ZERO, muda-o.
		else
			pRow->Update( node );
    }
	
	// Se nao encontrou a posicao do elemento e ele nao e' ZERO...
	else if ( ! IsZero( value ) )
    {
		node.col = col;
		pRow->Insert( node );
    }
	
	return( 1 );
}



/**************/
/*** GetVal ***/
//
//  Le um elemento da matriz, sem verificar fronteiras. O valor
// da linha (row) deve ser maior ou igual ao da coluna (col).
//
template<class TVar>
const TVar &
TPZSpMatrix<TVar>::GetVal(const int64_t row,const int64_t col ) const
{
	TPZLink<TPZNode> *pRow = &fElem[row];
	TPZNode        node;
	
	// Caso a lista esteja vazia, garante que (node.col != col).
	node.col = -1;
	
	// Procura pela posicao que esta ou deveria estar o elemento.
	pRow->Head();
	while ( pRow->Get( &node ) && (node.col < col) )
		pRow->Next();
	
	// Se encontrou a posicao do elemento...
	if ( node.col == col )
		return( pRow->GetNode()->elem );
	else {
        static TVar zero;
        zero = TVar(0);
		return( zero );
	}
}




/******** Operacoes com matrizes FULL  ********/

/******************/
/*** Operator = ***/

template<class TVar>
TPZSpMatrix<TVar> &
TPZSpMatrix<TVar>::operator=(const TPZSpMatrix<TVar> &A )
{
	delete [] fElem;
	fCopy( &A );
	return( *this );
}



/******************/
/*** Operator + ***/

template<class TVar>
TPZSpMatrix<TVar>
TPZSpMatrix<TVar>::operator+(const TPZSpMatrix<TVar> &A ) const
{
	TPZSpMatrix<TVar> res( *this );
	if ( ! res.fAdd( &A ) )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Operator+ (TPZSpMatrix&) <incompatible dimensions>" );
	
	return( res );
}



/******************/
/*** Operator - ***/

template<class TVar>
TPZSpMatrix<TVar>
TPZSpMatrix<TVar>::operator-(const TPZSpMatrix<TVar> &A ) const
{
	TPZSpMatrix<TVar> res( *this );
	if ( ! res.fSub( &A ) )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Operator+( TPZSpMatrix ) <incompatible dimensions>" );
	
	return( res );
}

/*******************/
/*** Operator += ***/

template<class TVar>
TPZSpMatrix<TVar> &
TPZSpMatrix<TVar>::operator+=(const TPZSpMatrix<TVar> &A )
{
	if ( ! fAdd( &A ) )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "operator+( TPZSpMatrix ) <incompatible dimensions>" );
	return( *this );
}



/*******************/
/*** Operator -= ***/

template<class TVar>
TPZSpMatrix<TVar> &
TPZSpMatrix<TVar>::operator-=(const TPZSpMatrix<TVar> &A )
{
	if ( ! fSub( &A ) )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "operator-( TPZSpMatrix ) <incompatible dimensions>" );
	return( *this );
}

/******** Operacoes com valores NUMERICOS ********/
/*****************************/
/*** Operator * ( REAL ) ***/
template<class TVar>
TPZSpMatrix<TVar>
TPZSpMatrix<TVar>::operator*(const TVar value ) const
{
	TPZSpMatrix<TVar> res( *this );
	res.fMult( value );
	return( res );
}

/******************************/
/*** Operator *= ( REAL ) ***/

template<class TVar>
TPZSpMatrix<TVar> &
TPZSpMatrix<TVar>::operator*=(const TVar value )
{
	if ( IsZero( value ) )
		return( Reset() );
	
	fMult( value );
	return( *this );
}



/******** Metodos genericos ********/

/*************/
/*** Reset ***/
template<class TVar>
TPZSpMatrix<TVar> &
TPZSpMatrix<TVar>::Reset()
{
	for ( int64_t i = 0; i < this->Rows(); i++ )
		fElem[i].Clear();
	return( *this );
}

/**************/
/*** Resize ***/
//
// Muda as dimensoes da matriz, mas matem seus valores antigos. Novas
// posicoes sao criadas com ZEROS.
//
template<class TVar>
int
TPZSpMatrix<TVar>::Resize(const int64_t newRows,const int64_t newCols )
{
	if ( newRows == this->Rows() )
		return( 1 );
	
	// Cria nova matrix.
	TPZLink<TPZNode> *newDiag = new TPZLink<TPZNode>[ newRows ] ;
	
	// Copia os elementos para a nova matriz.
	int64_t min = MIN( newRows, this->Rows() );
	for ( int64_t i = 0; i < min; i++ )
		newDiag[i] = fElem[i];
	
	// Descarta a matriz antiga e valida a nova matriz.
	delete( fElem );
	fElem = newDiag;
	// fRow  = fCol = newRows;
	this->fRow=newRows;this->fCol=newCols;
	return( 1 );
}

/*************/
/*** Redim ***/
//
// Muda as dimensoes da matriz e ZERA seus elementos.
//
template<class TVar>
int
TPZSpMatrix<TVar>::Redim(const int64_t newRows,const int64_t newCols )
{
	this->fCol = newCols;
	delete [] fElem;
	fElem = new TPZLink<TPZNode>[ newRows ];
	//fRow  = fCol = newRows;
	this->fRow=newRows;
	return( 1 );
}

/*************/
/*** Zero ***/
template<class TVar>
int
TPZSpMatrix<TVar>::Zero()
{
	
	delete [] fElem;
	fElem = new TPZLink<TPZNode>[ this->fRow ];
	this->fDecomposed = 0;
	return( 1 );
}

/************************** PROTECTED **************************/

/*****************/
/*** fAdd (mat) ***/

template<class TVar>
int
TPZSpMatrix<TVar>::fAdd(const TPZSpMatrix<TVar> *const A )
{
	if ( (this->Rows() != A->Rows()) || (this->Cols() != A->Cols()) )
		return( 0 );
	
	TPZNode mNode, aNode;
	TPZLink<TPZNode> *pm = &fElem[0];
	TPZLink<TPZNode> *pa = &A->fElem[0];
	
	for ( int64_t row = 0; row < this->Rows(); row++, pm++, pa++ )
    {
		// Soma uma linha.
		pm->Head();
		pa->Head();
		int64_t mOk = pm->Get( &mNode );
		int64_t aOk = pa->Get( &aNode );
		
		// Enquanto as duas linhas tiverem elementos...
		while ( mOk && aOk )
		{
			// Se as colunas forem iguais, soma os elementos.
			if ( mNode.col == aNode.col )
			{
				mNode.elem += aNode.elem;
				pm->Update( mNode );
				pm->Next();
				pa->Next();
				mOk = pm->Get( &mNode );
				aOk = pa->Get( &aNode );
			}
			
			// Se a coluna desta matriz for maior, insere o elemento
			//  da matriz A.
			else if ( mNode.col > aNode.col )
			{
				pm->Insert( aNode );
				
				// Volta a lista 'pm' `a posicao anterior.
				pm->Next();
				
				pa->Next();
				aOk = pa->Get( &aNode );
			}
			
			// Se a coluna da matriz A for maior, apenas avanca a
			//  lista 'pm' para a proxima posicao.
			else
			{
				pm->Next();
				mOk = pm->Get( &mNode );
			}
		}
		
		// Se apenas a linha da matriz A tiver elementos,
		// copia-a para esta matriz.
		if ( aOk )
		{
			while ( pa->Get( &aNode ) )
			{
				pm->Append( aNode );
				pa->Next();
			}
		}
    }
	
	return( 1 );
}

/*****************/
/*** fSub (mat) ***/

template<class TVar>
int
TPZSpMatrix<TVar>::fSub(const TPZSpMatrix<TVar> *const A )
{
	if ( (this->Rows() != A->Rows()) || (this->Cols() != A->Cols()) )
		return( 0 );
	
	TPZNode mNode, aNode;
	TPZLink<TPZNode> *pm = &fElem[0];
	TPZLink<TPZNode> *pa = &A->fElem[0];
	
	for ( int64_t row = 0; row < this->Rows(); row++, pm++, pa++ )
    {
		// Soma uma linha.
		pm->Head();
		pa->Head();
		int64_t mOk = pm->Get( &mNode );
		int64_t aOk = pa->Get( &aNode );
		
		// Enquanto as duas linhas tiverem elementos...
		while ( mOk && aOk )
		{
			// Se as colunas forem iguais, soma os elementos.
			if ( mNode.col == aNode.col )
			{
				mNode.elem -= aNode.elem;
				pm->Update( mNode );
				pm->Next();
				pa->Next();
				mOk = pm->Get( &mNode );
				aOk = pa->Get( &aNode );
			}
			
			// Se a coluna desta matriz for maior, insere o elemento
			//  com sinal trocado.
			else if ( mNode.col > aNode.col )
			{
				aNode.elem = -aNode.elem;
				pm->Insert( aNode );
				
				// Volta a lista 'pm' `a posicao anterior.
				pm->Next();
				
				pa->Next();
				aOk = pa->Get( &aNode );
			}
			
			// Se a coluna da matriz A for maior, apenas avanca a
			//  lista 'pm' para a proxima posicao.
			else
			{
				pm->Next();
				mOk = pm->Get( &mNode );
			}
		}
		
		// Se apenas a linha da matriz A tiver elementos,
		// copia-a para esta matriz (com sinal invertido).
		if ( aOk )
		{
			while ( pa->Get( &aNode ) )
			{
				aNode.elem = -aNode.elem;
				pm->Append( aNode );
				pa->Next();
			}
		}
    }
	
	return( 1 );
}

/*******************/
/*** fCopy (mat) ***/
//
//  Aloca as linhas e copia os elementos.
//
template<class TVar>
int
TPZSpMatrix<TVar>::fCopy(const TPZSpMatrix<TVar> *const A )
{
	this->fCol  = A->Cols();
	this->fRow  = A->Rows();
	fElem = new TPZLink<TPZNode>[ this->fRow ];
	
	if ( fElem == NULL )
		return( 0 );
	
	TPZLink<TPZNode> *pm = &fElem[0];
	TPZLink<TPZNode> *pa = &A->fElem[0];
	for ( int64_t i = 0; i < this->fRow; i++ )
		*pm++ = *pa++;
	
	return( 1 );
}

/*********************/
/*** fMult (value) ***/
template<class TVar>
int
TPZSpMatrix<TVar>::fMult(const TVar value )
{
	TPZNode        node;
	TPZLink<TPZNode> *pm = &fElem[0];
	
	for ( int64_t row = 0; row < this->Rows(); row++, pm++ )
    {
		pm->Head();
		while ( pm->Get( &node ) )
		{
			node.elem *= value;
			pm->Update( node );
			pm->Next();
		}
    }
	
	return( 1 );
}


/*************************** PRIVATE ***************************/


/****************/
/*** Prod Esc ***/
template<class TVar>
REAL
TPZSpMatrix<TVar>::ProdEsc( TPZLink<TPZNode> *row_i, TPZLink<TPZNode> *row_j,
					 int64_t k )
{
	TVar prod = 0.0;
	
	TPZNode node_i, node_j;
	
	row_i->Head();
	row_j->Head();
	
	int64_t again_i = row_i->Get( &node_i ) && (node_i.col < k);
	int64_t again_j = row_j->Get( &node_j ) && (node_j.col < k);
	
	while ( again_i && again_j )
    {
		if ( node_i.col > node_j.col )
		{
			row_j->Next();
			again_j = row_j->Get( &node_j ) && (node_j.col < k);
		}
		else if ( node_i.col < node_j.col )
		{
			row_i->Next();
			again_i = row_i->Get( &node_i ) && (node_i.col < k);
		}
		else if ( node_i.col == node_j.col )
		{
			prod += node_i.elem * node_j.elem;
			row_i->Next();
			row_j->Next();
			again_i = row_i->Get( &node_i ) && (node_i.col < k);
			again_j = row_j->Get( &node_j ) && (node_j.col < k);
		}
    }
	
	// Garante que a lista 'row_i' esteja na coluna 'k'.
	while ( again_i )
    {
		row_i->Next();
		again_i = row_i->Get( &node_i ) && (node_i.col < k);
    }
	
	// Garante que a lista 'row_j' esteja na coluna 'k'.
	while ( again_j )
    {
		row_j->Next();
		again_j = row_j->Get( &node_j ) && (node_j.col < k);
    }
	
	return( prod );
}

template<class TVar>
void TPZSpMatrix<TVar>::MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
						  const REAL alpha,const REAL beta,const int opt) const  {
	if ((!opt && this->Cols() != x.Rows()) || this->Rows() != x.Rows())
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "TPZSpMatrix::MultAdd <matrixs with incompatible dimensions>" );
	if(x.Cols() != y.Cols() || x.Cols() != z.Cols() || x.Rows() != y.Rows() || x.Rows() != z.Rows()) {
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"TPZSpMatrix::MultAdd incompatible dimensions\n");
	}
	int64_t rows = this->Rows();
	int64_t xcols = x.Cols();
	int64_t ic, r;
	this->PrepareZ(y,z,beta,opt);
	TVar val;
	for (ic = 0; ic < xcols; ic++) {
		if(!opt) {
			const REAL* firstel = &x.g(0,ic);
			TPZNode *currentnode;
			for ( r = 0; r < rows; r++ ) {
				TPZLink<TPZNode> *runner = &fElem[r];
				if(!runner->Head()) continue;
				val = 0.;
				currentnode = runner->GetNode();
				while(currentnode) {
					val += *(firstel+(currentnode->col)) * currentnode->elem;
					runner->Next();
					currentnode = runner->GetNode();
				}
				z(r,ic) += alpha*val;
			}
		} else {
			TVar * firstelz = &z(0,ic);
			TPZNode *currentnode;
			for (r = 0; r<rows; r++) {
				TVar elx = x.g(r,ic);
				TPZLink<TPZNode> *runner = &fElem[r];
				if(!runner->Head()) continue;
				do {
					currentnode = runner->GetNode();
					*(firstelz+(currentnode->col)) += alpha* elx * currentnode->elem;
				} while (runner->Next());
			}
		}
	}
}


#ifdef OOPARLIB

template<class TVar>
int TPZSpMatrix<TVar>::Unpack( TReceiveStorage *buf ){
	TMatrix::Unpack(buf);
	int64_t rows;
	buf->UpkInt(&rows);
	Redim(rows);
	int64_t nelem;
	int64_t col;
	TVar val;
	for(int64_t i=0;i<rows;i++) {
		buf->UpkInt(&nelem);
		buf->UpkDouble(&val);
		buf->UpkInt(&col);
		PutVal(i,col,val);
	}
	return 1;
}


template<class TVar>
TSaveable *TPZSpMatrix<TVar>::CreateInstance(TReceiveStorage *buf) {
	TPZSpMatrix<TVar> *m = new TPZSpMatrix<TVar>();
	m->Unpack(buf);
	return m;
}

template<class TVar>
int TPZSpMatrix<TVar>::Pack( TSendStorage *buf ) const {
	TMatrix::Pack(buf);
	TPZNode        node;
	TPZLink<TPZNode> *pm = &fElem[0];
	int64_t rows = Rows();
	buf->PkInt(&rows);
	for ( int64_t row = 0; row < rows; row++, pm++ )
    {
		int64_t numel = 0;
		pm->Head();
		while ( pm->Get( &node ) )
		{
			numel++;
			pm->Next();
		}
		buf->PkInt(&numel);
		pm->Head();
		while ( pm->Get( &node ) )
		{
			// guardar os dados de node
			buf->PkDouble(&node.elem);
			buf->PkInt(&node.col);
			pm->Next();
		}
    }
	return 1;
}

template<class TVar>
int TPZSpMatrix<TVar>::DerivedFrom(const int64_t Classid) const {
	return TMatrix::DerivedFrom(Classid);
}

template<class TVar>
int TPZSpMatrix<TVar>::DerivedFrom(const char *classname) const {
	
	if(!strcmp(ClassName(),classname)) return 1;
	return TMatrix::DerivedFrom(classname);
}

#endif

template<class TVar>
int TPZSpMatrix<TVar>::ClassId() const{
    return Hash("TPZSpMatrix") ^ TPZMatrix<TVar>::ClassId() << 1;
}
template class TPZSpMatrix<REAL>;
