/**
 * @file
 * @brief Contains the implementation of the TPZSSpMatrix methods.
 */
//
//
// Author: MISAEL LUIS SANTANA MANDUJANO.
//
// File:   tsespmat.cc
//
// Class:  TPZSSpMatrix
//
// Obs.:   Implementa matrizes esparsas simetricas. A implementacao
//         e' feita atraves de listas ligadas.
//
// Versao: 04 / 1996.
//

#include <math.h>
#include <stdlib.h>

#include "pzfmatrix.h"
#include "pzsespmat.h" 
//#include "pzerror.h"

#include <sstream>
#include "pzlog.h"
#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.matrix.tpzsspmatrix"));
#endif

using namespace std;

/*******************/
/*** TPZSSpMatrix ***/

/**************************** PUBLIC ****************************/

/*******************/
/*** Constructor ***/

template<class TVar>
TPZSSpMatrix<TVar>::TPZSSpMatrix(const TPZSSpMatrix<TVar> &A )
: TPZMatrix<TVar>( A ), fMat( A.fMat )
{
	
    this->fDecomposed  = A.fDecomposed;
    this->fDefPositive = A.fDefPositive;
}

/******** Operacoes com matrizes ESPARSAS SIMETRICAS  ********/

/******************/
/*** Operator = ***/

template<class TVar>
TPZSSpMatrix<TVar> &
TPZSSpMatrix<TVar>::operator=(const TPZSSpMatrix<TVar> &A )
{
    fMat = A.fMat;
    this->fRow = this->fCol = A.Dim();
    this->fDecomposed  = A.fDecomposed;
    this->fDefPositive = A.fDefPositive;
    return( *this );
}



/******************/
/*** Operator + ***/

template<class TVar>
TPZSSpMatrix<TVar>
TPZSSpMatrix<TVar>::operator+(const TPZSSpMatrix<TVar> &A ) const
{
    if ( this->Dim() != A.Dim() )
        TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Operator+ (TPZSSpMatrix&) <incompatible dimensions>" );
    TPZSSpMatrix res( *this );
    res.fMat.fAdd( &A.fMat );
    return( res );
}



/******************/
/*** Operator - ***/

template<class TVar>
TPZSSpMatrix<TVar>
TPZSSpMatrix<TVar>::operator-(const TPZSSpMatrix<TVar> &A ) const
{
    if ( this->Dim() != A.Dim() )
        TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Operator-( TPZSSpMatrix& ) <incompatible dimensions>" );
    
    TPZSSpMatrix res( *this );
    res.fMat.fSub( &A.fMat );
    return( res );
}



/*******************/
/*** Operator += ***/

template<class TVar>
TPZSSpMatrix<TVar> &
TPZSSpMatrix<TVar>::operator+=(const TPZSSpMatrix<TVar> &A )
{
    if ( this->Dim() != A.Dim() )
        TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "operator+( TPZSSpMatrix ) <incompatible dimensions>" );
    
    fMat.fAdd( &A.fMat );
    this->fDecomposed = 0;
    return( *this );
}



/*******************/
/*** Operator -= ***/

template<class TVar>
TPZSSpMatrix<TVar> &
TPZSSpMatrix<TVar>::operator-=(const TPZSSpMatrix<TVar> &A )
{
    if ( this->Dim() != A.Dim() )
        TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "operator+( TPZSSpMatrix ) <incompatible dimensions>" );
    
    fMat.fSub( &A.fMat );
    this->fDecomposed = 0;
    return( *this );
}



/******** Operacoes com matrizes GENERICAS ********/

// Estas operacoes com matrizes genericas, usam a parte triangular
// inferior da maior matriz quadrada de A. Ex.:
//
//  Se A = 01 02 03 04   A matriz usada sera':  01 05 09
//         05 06 07 08                          05 06 10
//         09 10 11 12                          09 10 11
//



/******** Operacoes com valores NUMERICOS ********/



/*****************************/
/*** Operator * ( REAL ) ***/

template<class TVar>
TPZSSpMatrix<TVar>
TPZSSpMatrix<TVar>::operator*(const TVar value ) const
{
    TPZSSpMatrix<TVar> res( *this );
    res.fMat.fMult( value );
    return( res );
}




/******************************/
/*** Operator *= ( TVar ) ***/

template<class TVar>
TPZSSpMatrix<TVar> &
TPZSSpMatrix<TVar>::operator*=(const TVar value )
{
    fMat.fMult( value );
    this->fDecomposed = 0;
    return( *this );
}



/******** Metodos genericos ********/

/**************************/
/*** Decompose Cholesky ***/
template<class TVar>
int
TPZSSpMatrix<TVar>::Decompose_Cholesky(std::list<int> &singular)
{
    return Decompose_Cholesky();
}

template<class TVar>
int
TPZSSpMatrix<TVar>::Decompose_Cholesky()
{
    
    if ( this->fDecomposed && this->fDecomposed != ECholesky)  TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Decompose_Cholesky <Matrix already Decomposed with different algorithm>" );
    else if(this->fDecomposed) return 0;
    
    
    TPZSpMatrix<REAL>::TPZNode node_k;
    TPZLink<TPZSpMatrix<REAL>::TPZNode> *row_end = &fMat.fElem[ this->Dim() ];
    TPZLink<TPZSpMatrix<REAL>::TPZNode> *row_k   = &fMat.fElem[ 0 ];
    
    for ( int k = 0; k < this->Dim(); k++, row_k++ )
    {
        // Faz sum = SOMA( A(k,p) * A(k,p) ), p = 1, ..., k-1.
        //
        REAL sum = 0.0;
        row_k->Head();
        node_k.col = -1;
        for (; row_k->Get( &node_k ) && (node_k.col < k);
             row_k->Next() )
            sum += node_k.elem * node_k.elem;
        
        // Faz A(k,k) = sqrt( A(k,k) - sum ).
        //
        node_k.elem -= sum;
        if ( (node_k.col != k) || (node_k.elem < 1.e-10) )
            return( 0 );
        //	TPZMatrix::Error(__PRETTY_FUNCTION__, "Cholesky_Decomposition ",
        //	       "<matrix is not positive definite>" );
        if ( node_k.elem < 1.e-10 )
            return( 0 );
        //	TPZMatrix::Error(__PRETTY_FUNCTION__, "Cholesky_Decomposition <sqrt of negative number>");
        TVar pivot = node_k.elem = sqrt( node_k.elem );
        row_k->Update( node_k );
        
        // Loop para i = k+1 ... Dim().
        //
        TPZSpMatrix<REAL>::TPZNode node_i;
        TPZLink<TPZSpMatrix<REAL>::TPZNode> *row_i = row_k + 1;
        for ( ; row_i != row_end; row_i++ )
        {
            // Faz sum = SOMA( A(i,p) * A(k,p) ), p = 1, ..., k-1.
            //
            sum = ProdEsc( row_k, row_i, k );
            
            // Faz A(i,k) = (A(i,k) - sum) / A(k,k)
            //
            node_i.col = -1;
            row_i->Get( &node_i );
            if ( node_i.col != k )
            {
                node_i.elem = -sum / pivot;
                node_i.col = k;
                row_i->Insert( node_i );
            }
            else
            {
                node_i.elem = (node_i.elem - sum) / pivot;
                row_i->Update( node_i );
            }
        }
    }
    
    this->fDecomposed  = 1;
    this->fDefPositive = 1;
    return( 1 );
}



/**********************/
/*** Decompose LDLt ***/
template<class TVar>
int
TPZSSpMatrix<TVar>::Decompose_LDLt(std::list<int> &singular)
{
    return Decompose_LDLt();
}

template<class TVar>
int
TPZSSpMatrix<TVar>::Decompose_LDLt()
{
    if (  this->fDecomposed && this->fDecomposed != ELDLt)  TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Decompose_LDLt <Matrix already Decomposed with different method>" );
    else if(this->fDecomposed) return 0;
    
	
    TPZSpMatrix<REAL>::TPZNode node_k;
    TPZLink<TPZSpMatrix<REAL>::TPZNode> row_aux;
    TPZLink<TPZSpMatrix<REAL>::TPZNode> *row_k   = &fMat.fElem[ 0 ];
    TPZLink<TPZSpMatrix<REAL>::TPZNode> *row_end = &fMat.fElem[ this->Dim() ];
    
    for ( int k = 0; k < this->Dim(); k++, row_k++ )
    {
        // Faz aux(p) = A(p,p) * A(k,p), p = 1, ..., k-1.
        //
        row_k->Head();
        row_aux.Clear();
        for ( ; row_k->Get( &node_k ) && (node_k.col < k);
             row_k->Next() )
        {
            node_k.elem *= GetVal( node_k.col, node_k.col );
            if ( node_k.elem!=REAL(0.0) )
                row_aux.Append( node_k );
        }
        
        // Faz sum = SOMA( A(k,p) * aux(p) ), p = 1, ..., k-1.
        //
        TVar sum = ProdEsc( row_k, &row_aux, k );
        
        
        // Faz A(k,k) = A(k,k) - sum;
        //
        node_k.col = -1;
        row_k->Get( &node_k );
        if ( node_k.col != k )
        {
            // O elemento antigo A(k,k) e' ZERO...
            node_k.col  = k;
            node_k.elem = -sum;
            if ( IsZero( node_k.elem ) )
                TPZMatrix<REAL> ::Error(__PRETTY_FUNCTION__, "LDLt_Decomposition <diag element is zero>" );
            row_k->Insert( node_k );
        }
        else
        {
            // O elemento antigo A(k,k) nao e' ZERO...
            node_k.elem -= sum;
            if ( IsZero( node_k.elem ) )
                TPZMatrix<REAL> ::Error(__PRETTY_FUNCTION__, "LDLt_Decomposition <diag element is zero>" );
            row_k->Update( node_k );
        }
        TVar pivot = node_k.elem;
        
        // Loop para i = k+1, ..., Dim().
        //
        TPZSpMatrix<REAL>::TPZNode node_i;
        TPZLink<TPZSpMatrix<REAL>::TPZNode> *row_i = row_k + 1;
        for ( ; row_i != row_end; row_i++ )
        {
            // Faz sum = SOMA( A(i,p) * aux(p) ), p = 1, ..., k-1.
            //
            sum = ProdEsc( row_i, &row_aux, k );
            
            // Faz A(i,k) = (A(i,k) - sum) / A(k,k).
            //
            node_i.col = -1;
            row_i->Get( &node_i );
            if ( node_i.col != k )
            {
                node_i.col = k;
                node_i.elem = -sum / pivot;
                row_i->Insert( node_i );
            }
            else
            {
                node_i.elem = (node_i.elem - sum) / pivot;
                row_i->Update( node_i );
            }
        }
    }
    
    this->fDecomposed  = 1;
    this->fDefPositive = 0;
    return( 1 );
}



/*********************/
/*** Subst Forward ***/
//
//  Faz Ax = b, onde A e' triangular inferior.
//
template<class TVar>
int
TPZSSpMatrix<TVar>::Subst_Forward( TPZFMatrix<TVar> *B ) const
{
    if ( (B->Rows() != this->Dim()) || !this->fDecomposed )
        return( 0 );
    
    TPZLink<TPZSpMatrix<REAL>::TPZNode> *row_k = &fMat.fElem[0];
    for ( int k = 0; k < this->Dim(); k++, row_k++ )
        for ( int j = 0; j < B->Cols(); j++ )
        {
            // Faz sum = SOMA( A[k,i] * B[i,j] ), para i = 1, ..., k-1.
            //
            TVar sum = 0.0;
            TPZSpMatrix<REAL>::TPZNode  ni;
            row_k->Head();
            for ( ; row_k->Get( &ni ) && (ni.col < k ); row_k->Next() )
                sum += ni.elem * B->GetVal( ni.col, j );
            
            // E' assumido que A[k,k] nunca e' ZERO, pois a matriz
            //  ja' foi decomposta.
            // A[k,k] == ni.elem, se ni.col == k.
            
            // Faz B[k,j] = (B[k,j] - sum) / A[k,k].
            //
            B->PutVal( k, j, (B->GetVal( k, j ) - sum) / ni.elem );
        }
    
    return( 1 );
}



/***********************/
/*** Subst L Forward ***/
//
//  Faz a "Forward substitution" assumindo que os elementos
//   da diagonal sao todos iguais a 1.
//

template<class TVar>
int
TPZSSpMatrix<TVar>::Subst_LForward( TPZFMatrix<TVar> *B ) const
{
    if ( (B->Rows() != this->Dim()) || !this->fDecomposed )
        return( 0 );
    
    TPZLink<TPZSpMatrix<REAL>::TPZNode> *row_k = &fMat.fElem[0];
    for ( int k = 0; k < this->Dim(); k++, row_k++ )
        for ( int j = 0; j < B->Cols(); j++ )
        {
            // Faz sum = SOMA( A[k,i] * B[i,j] ), para i = 1, ..., k-1.
            //
            TVar sum = 0.0;
            TPZSpMatrix<REAL>::TPZNode  ni;
            row_k->Head();
            for ( ; row_k->Get( &ni ) && (ni.col < k ); row_k->Next() )
                sum += ni.elem * B->GetVal( ni.col, j );
            
            // Faz b[k] = (b[k] - sum).
            //
            B->PutVal( k, j, B->GetVal( k, j ) - sum );
        }
    
    return( 1 );
}



/******************/
/*** Subst Diag ***/
//
//  Faz Ax = B, sendo que A e' assumida ser uma matriz diagonal.
//
template<class TVar>
int
TPZSSpMatrix<TVar>::Subst_Diag( TPZFMatrix<TVar> *B ) const
{
    // So' e' permitido fazer esta substituicao se a matriz foi
    //  decomposta, pois desta forma se garante que os elementos
    //  da diagonal sao diferentes de ZERO. O que elimina uma
    //  verificacao.
    //
    if ( (B->Rows() != this->Dim()) || !this->fDecomposed )
        return( 0 );
    
    TPZSpMatrix<REAL>::TPZNode node_kk;
    TPZLink<TPZSpMatrix<REAL>::TPZNode> *row_k = &fMat.fElem[0];
    for ( int k = 0; k < this->Dim(); k++, row_k++ )
    {
        row_k->GetLast( &node_kk );
        for ( int j = 0; j < B->Cols(); j++ )
            B->PutVal( k, j, B->GetVal( k, j ) / node_kk.elem );
    }
    return( 1 );
}




/**************************** PRIVATE ****************************/



/*************/
/*** Error ***/
/*int
 TPZSSpMatrix::Error(const char *msg1,const char *msg2 ) 
 {
 ostringstream out;
 out << "TPZSSpMatrix::" << msg1 << msg2 << ".\n";
 //pzerror.show();
 LOGPZ_ERROR (logger, out.str().c_str());
 DebugStop();
 return 0;
 }*/



/****************/
/*** Prod Esc ***/
template<class TVar>
TVar
TPZSSpMatrix<TVar>::ProdEsc( TPZLink<TPZSpMatrix<REAL>::TPZNode> *row_i,
                      TPZLink<TPZSpMatrix<REAL>::TPZNode> *row_j, int k )
{
    TVar prod = 0.0;
    
    TPZSpMatrix<REAL>::TPZNode node_i, node_j;
    
    row_i->Head();
    row_j->Head();
    
    int again_i = row_i->Get( &node_i ) && (node_i.col < k);
    int again_j = row_j->Get( &node_j ) && (node_j.col < k);
    
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

#ifdef OOPARLIB

template<class TVar>
int TPZSSpMatrix<TVar>::Unpack (TReceiveStorage *buf ){
    TSimMatrix::Unpack(buf);
    long class_id;
    buf->UpkLong(&class_id);
    if(!fMat.DerivedFrom(class_id)) DebugStop();
    fMat.Unpack(buf);
    return 1;
}


template<class TVar>
TSaveable *TPZSSpMatrix<TVar>::Restore(TReceiveStorage *buf) {
    TPZSSpMatrix<TVar> *m = new TPZSSpMatrix<TVar>();
    m->Unpack(buf);
    return m;
}

template<class TVar>
int TPZSSpMatrix<TVar>::Pack( TSendStorage *buf ) const {
    TSimMatrix::Pack(buf);
    fMat.Pack(buf);
    return 1;
}

template<class TVar>
int TPZSSpMatrix<TVar>::DerivedFrom(const long Classid) const {
    if(Classid == GetClassID()) return 1;
    return TSimMatrix::DerivedFrom(Classid);
}

template<class TVar>
int TPZSSpMatrix<TVar>::DerivedFrom(const char *classname) const {
    if(!strcmp(ClassName(),classname)) return 1;
    return TSimMatrix::DerivedFrom(classname);
}


#endif

