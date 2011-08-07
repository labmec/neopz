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

TPZSSpMatrix::TPZSSpMatrix(const TPZSSpMatrix &A )
: TPZMatrix( A ), fMat( A.fMat )
{
    fDecomposed  = A.fDecomposed;
    fDefPositive = A.fDefPositive;
}

/******** Operacoes com matrizes ESPARSAS SIMETRICAS  ********/

/******************/
/*** Operator = ***/

TPZSSpMatrix &
TPZSSpMatrix::operator=(const TPZSSpMatrix &A )
{
    fMat = A.fMat;
    fRow = fCol = A.Dim();
    fDecomposed  = A.fDecomposed;
    fDefPositive = A.fDefPositive;
    return( *this );
}



/******************/
/*** Operator + ***/

TPZSSpMatrix
TPZSSpMatrix::operator+(const TPZSSpMatrix &A ) const
{
    if ( Dim() != A.Dim() )
        TPZMatrix::Error(__PRETTY_FUNCTION__, "Operator+ (TPZSSpMatrix&) <incompatible dimensions>" );
    TPZSSpMatrix res( *this );
    res.fMat.fAdd( &A.fMat );
    return( res );
}



/******************/
/*** Operator - ***/

TPZSSpMatrix
TPZSSpMatrix::operator-(const TPZSSpMatrix &A ) const
{
    if ( Dim() != A.Dim() )
        TPZMatrix::Error(__PRETTY_FUNCTION__, "Operator-( TPZSSpMatrix& ) <incompatible dimensions>" );
    
    TPZSSpMatrix res( *this );
    res.fMat.fSub( &A.fMat );
    return( res );
}



/*******************/
/*** Operator += ***/

TPZSSpMatrix &
TPZSSpMatrix::operator+=(const TPZSSpMatrix &A )
{
    if ( Dim() != A.Dim() )
        TPZMatrix::Error(__PRETTY_FUNCTION__, "operator+( TPZSSpMatrix ) <incompatible dimensions>" );
    
    fMat.fAdd( &A.fMat );
    fDecomposed = 0;
    return( *this );
}



/*******************/
/*** Operator -= ***/

TPZSSpMatrix &
TPZSSpMatrix::operator-=(const TPZSSpMatrix &A )
{
    if ( Dim() != A.Dim() )
        TPZMatrix::Error(__PRETTY_FUNCTION__, "operator+( TPZSSpMatrix ) <incompatible dimensions>" );
    
    fMat.fSub( &A.fMat );
    fDecomposed = 0;
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

TPZSSpMatrix
TPZSSpMatrix::operator*(const REAL value ) const
{
    TPZSSpMatrix res( *this );
    res.fMat.fMult( value );
    return( res );
}




/******************************/
/*** Operator *= ( REAL ) ***/

TPZSSpMatrix &
TPZSSpMatrix::operator*=(const REAL value )
{
    fMat.fMult( value );
    fDecomposed = 0;
    return( *this );
}



/******** Metodos genericos ********/

/**************************/
/*** Decompose Cholesky ***/
int
TPZSSpMatrix::Decompose_Cholesky(std::list<int> &singular)
{
    return Decompose_Cholesky();
}

int
TPZSSpMatrix::Decompose_Cholesky()
{
    
    if ( fDecomposed && fDecomposed != ECholesky)  TPZMatrix::Error(__PRETTY_FUNCTION__, "Decompose_Cholesky <Matrix already Decomposed with different algorithm>" );
    else if(fDecomposed) return 0;
    
    
    TPZSpMatrix::TPZNode node_k;
    TPZLink<TPZSpMatrix::TPZNode> *row_end = &fMat.fElem[ Dim() ];
    TPZLink<TPZSpMatrix::TPZNode> *row_k   = &fMat.fElem[ 0 ];
    
    for ( int k = 0; k < Dim(); k++, row_k++ )
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
        REAL pivot = node_k.elem = sqrt( node_k.elem );
        row_k->Update( node_k );
        
        // Loop para i = k+1 ... Dim().
        //
        TPZSpMatrix::TPZNode node_i;
        TPZLink<TPZSpMatrix::TPZNode> *row_i = row_k + 1;
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
    
    fDecomposed  = 1;
    fDefPositive = 1;
    return( 1 );
}



/**********************/
/*** Decompose LDLt ***/
int
TPZSSpMatrix::Decompose_LDLt(std::list<int> &singular)
{
    return Decompose_LDLt();
}

int
TPZSSpMatrix::Decompose_LDLt()
{
    if (  fDecomposed && fDecomposed != ELDLt)  TPZMatrix::Error(__PRETTY_FUNCTION__, "Decompose_LDLt <Matrix already Decomposed with different method>" );
    else if(fDecomposed) return 0;
    
    TPZSpMatrix::TPZNode node_k;
    TPZLink<TPZSpMatrix::TPZNode> row_aux;
    TPZLink<TPZSpMatrix::TPZNode> *row_k   = &fMat.fElem[ 0 ];
    TPZLink<TPZSpMatrix::TPZNode> *row_end = &fMat.fElem[ Dim() ];
    
    for ( int k = 0; k < Dim(); k++, row_k++ )
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
        REAL sum = ProdEsc( row_k, &row_aux, k );
        
        
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
                TPZMatrix::Error(__PRETTY_FUNCTION__, "LDLt_Decomposition <diag element is zero>" );
            row_k->Insert( node_k );
        }
        else
        {
            // O elemento antigo A(k,k) nao e' ZERO...
            node_k.elem -= sum;
            if ( IsZero( node_k.elem ) )
                TPZMatrix::Error(__PRETTY_FUNCTION__, "LDLt_Decomposition <diag element is zero>" );
            row_k->Update( node_k );
        }
        REAL pivot = node_k.elem;
        
        // Loop para i = k+1, ..., Dim().
        //
        TPZSpMatrix::TPZNode node_i;
        TPZLink<TPZSpMatrix::TPZNode> *row_i = row_k + 1;
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
    
    fDecomposed  = 1;
    fDefPositive = 0;
    return( 1 );
}



/*********************/
/*** Subst Forward ***/
//
//  Faz Ax = b, onde A e' triangular inferior.
//
int
TPZSSpMatrix::Subst_Forward( TPZFMatrix *B ) const
{
    if ( (B->Rows() != Dim()) || !fDecomposed )
        return( 0 );
    
    TPZLink<TPZSpMatrix::TPZNode> *row_k = &fMat.fElem[0];
    for ( int k = 0; k < Dim(); k++, row_k++ )
        for ( int j = 0; j < B->Cols(); j++ )
        {
            // Faz sum = SOMA( A[k,i] * B[i,j] ), para i = 1, ..., k-1.
            //
            REAL sum = 0.0;
            TPZSpMatrix::TPZNode  ni;
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
int
TPZSSpMatrix::Subst_LForward( TPZFMatrix *B ) const
{
    if ( (B->Rows() != Dim()) || !fDecomposed )
        return( 0 );
    
    TPZLink<TPZSpMatrix::TPZNode> *row_k = &fMat.fElem[0];
    for ( int k = 0; k < Dim(); k++, row_k++ )
        for ( int j = 0; j < B->Cols(); j++ )
        {
            // Faz sum = SOMA( A[k,i] * B[i,j] ), para i = 1, ..., k-1.
            //
            REAL sum = 0.0;
            TPZSpMatrix::TPZNode  ni;
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
int
TPZSSpMatrix::Subst_Diag( TPZFMatrix *B ) const
{
    // So' e' permitido fazer esta substituicao se a matriz foi
    //  decomposta, pois desta forma se garante que os elementos
    //  da diagonal sao diferentes de ZERO. O que elimina uma
    //  verificacao.
    //
    if ( (B->Rows() != Dim()) || !fDecomposed )
        return( 0 );
    
    TPZSpMatrix::TPZNode node_kk;
    TPZLink<TPZSpMatrix::TPZNode> *row_k = &fMat.fElem[0];
    for ( int k = 0; k < Dim(); k++, row_k++ )
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

REAL
TPZSSpMatrix::ProdEsc( TPZLink<TPZSpMatrix::TPZNode> *row_i,
                      TPZLink<TPZSpMatrix::TPZNode> *row_j, int k )
{
    REAL prod = 0.0;
    
    TPZSpMatrix::TPZNode node_i, node_j;
    
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

int TPZSSpMatrix::Unpack (TReceiveStorage *buf ){
    TSimMatrix::Unpack(buf);
    long class_id;
    buf->UpkLong(&class_id);
    if(!fMat.DerivedFrom(class_id)) DebugStop();
    fMat.Unpack(buf);
    return 1;
}



TSaveable *TPZSSpMatrix::Restore(TReceiveStorage *buf) {
    TPZSSpMatrix *m = new TPZSSpMatrix();
    m->Unpack(buf);
    return m;
}

int TPZSSpMatrix::Pack( TSendStorage *buf ) const {
    TSimMatrix::Pack(buf);
    fMat.Pack(buf);
    return 1;
}


int TPZSSpMatrix::DerivedFrom(const long Classid) const {
    if(Classid == GetClassID()) return 1;
    return TSimMatrix::DerivedFrom(Classid);
}

int TPZSSpMatrix::DerivedFrom(const char *classname) const {
    if(!strcmp(ClassName(),classname)) return 1;
    return TSimMatrix::DerivedFrom(classname);
}


#endif

