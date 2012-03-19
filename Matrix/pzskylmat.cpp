/**
 * @file
 * @brief Contains the implementation of the TPZSkylMatrix methods.
 */
//
// Author: MISAEL LUIS SANTANA MANDUJANO.
//
// File:   tskylmat.cc
//
// Class:  TPZSkylMatrix
//
// Obs.:   ESTA CLASSE GERENCIA MATRIZES DO TIPO SKYLINE.
//
// Versao: 4 / 1996.
//

#include <math.h>
#include <stdlib.h>

#ifdef BLAS
extern "C" {
#include <cblas.h>
}
#endif

/**
 * Commented out by Longhin
 * Compilation problem under MaCOSX OS.
 * Used by method TestSpeed which was also commented out.
 */
//#include <sys/timeb.h>

#include "pzfmatrix.h"
#include "pzskylmat.h"
//#include "pzerror.h"

#include <sstream>
#include "pzlog.h"
#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.matrix.tpzskylmatrix"));
#endif

using namespace std;

/*******************/
/*** TPZSkylMatrix ***/

/**************************** PUBLIC ****************************/

/*****************************/
/*** Construtor (int) ***/

template<class TVar>
TPZSkylMatrix<TVar>::TPZSkylMatrix(const int dim )
: TPZMatrix<TVar>( dim, dim ), fElem(dim+1), fStorage(0)
{
	
	// Inicializa a diagonal (vazia).
	fElem.Fill(0);
}
template<class TVar>
TPZSkylMatrix<TVar>::TPZSkylMatrix(const int dim, const TPZVec<int> &skyline )
: TPZMatrix<TVar>( dim, dim ), fElem(dim+1), fStorage(0)
{
	
	// Inicializa a diagonal (vazia).
	fElem.Fill(0);
	InitializeElem(skyline,fStorage,fElem);
}

template<class TVar>
void TPZSkylMatrix<TVar>::AddSameStruct(TPZSkylMatrix<TVar> &B, double k){
#ifdef DEBUG
	{
		int size = this->fElem.NElements();
		if(size != B.fElem.NElements()){
			PZError << "\nFATAL ERROR at " << __PRETTY_FUNCTION__ << "\n";
			PZError.flush();
			DebugStop();
		}
		for(int i = 0; i < size; i++){
			if((this->fElem[i]-this->fElem[0]) != (B.fElem[i]-B.fElem[0])){
				PZError << "\nFATAL ERROR at " << __PRETTY_FUNCTION__ << "\n";
				PZError.flush();
				DebugStop();
			}
		}
	}
#endif
	
	const int n = this->fStorage.NElements();
	for(int i = 0; i < n; i++) this->fStorage[i] += k*B.fStorage[i];
	
}

/**
 * @brief Updates the values of the matrix based on the values of the matrix
 */
template<class TVar>
void TPZSkylMatrix<TVar>::UpdateFrom(TPZAutoPointer<TPZMatrix<TVar> > mat)
{
    TPZMatrix<TVar> *matrix = mat.operator->();
    TPZSkylMatrix<TVar> *skylmat = dynamic_cast<TPZSkylMatrix<TVar> *>(matrix);
    if (!skylmat) {
        DebugStop();
    }
    if (fStorage.NElements() != skylmat->fStorage.NElements()) {
        DebugStop();
    }
    memcpy(&fStorage[0], &(skylmat->fStorage[0]) , fStorage.NElements()*sizeof(REAL));
    this->fDecomposed = skylmat->fDecomposed;
    this->fDefPositive = skylmat->fDefPositive;
}


template<class TVar>
void TPZSkylMatrix<TVar>::SetSkyline(const TPZVec<int> &skyline)
{
	fElem.Fill(0);
	InitializeElem(skyline,fStorage,fElem);
}
template<class TVar>
int TPZSkylMatrix<TVar>::NumElements(const TPZVec<int> &skyline) {
	int dim = skyline.NElements();
	int i,nelem=0;
	for(i=0; i<dim; i++) nelem += i-skyline[i]+1;
	return nelem;
}

template<class TVar>
void TPZSkylMatrix<TVar>::InitializeElem(const TPZVec<int> &skyline, TPZManVector<REAL> &storage, TPZVec<REAL *> &point) {
	int dim = skyline.NElements();
	int nel = NumElements(skyline);
	storage.Resize(nel);
	storage.Fill(0.);
	int i;
	point.Resize(dim+1);
	if(dim) {
		point[0] = &storage[0];
		point[dim] = &storage[0]+nel;
	} else {
		point[0] = 0;
	}
	for(i=1; i<dim+1; i++) point[i] = point[i-1]+(i-1)-skyline[i-1]+1;
}

/**
 Computes the highest skyline of both objects
 */
template<class TVar>
void TPZSkylMatrix<TVar>::ComputeMaxSkyline(const TPZSkylMatrix<TVar> &first, const TPZSkylMatrix<TVar> &second, TPZVec<int> &res) {
	
	if (first.Rows() != second.Rows()) {
		cout<<"ComputeMaxSkyline : incompatible dimension";
		return;
	}
	int i, dim = first.Rows();
	res.Resize(dim+1);
	
	for(i=1; i<dim+1; i++) {
		
		int aux = ( first.Size(i) > second.Size(i) ) ? first.Size(i) : second.Size(i);
		res[i] = i-aux-1;
	}
}

template<class TVar>
TVar &
TPZSkylMatrix<TVar>::operator()(const int r, const int c) {
	int row(r),col(c);
	if ( row > col ) this->Swap( &row, &col );
	
	// Indice do vetor coluna.
	int index = col - row;
	if ( index >= Size(col) ) {
		//Error("TPZSkylMatrix::operator()","Index out of range");
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"Index out of range");
		DebugStop();
	}
	return fElem[col][index];
}

template<class TVar>
TVar &TPZSkylMatrix<TVar>::s(const int row, const int col) {
	return operator()(row,col);
}

template<class TVar>
TVar &
TPZSkylMatrix<TVar>::operator()(const int r) {
	return operator()(r,r);
}



/**************/
/*** PutVal ***/
template<class TVar>
int
TPZSkylMatrix<TVar>::PutVal(const int r,const int c,const TVar & value )
{
	// inicializando row e col para trabalhar com a triangular superior
	int row(r),col(c);
	if ( row > col )
		this->Swap( &row, &col );
	
	// Indice do vetor coluna.
	int index = col - row;
	// Se precisar redimensionar o vetor.
	if ( index >= Size(col) && !IsZero(value)) {
		cout << "TPZSkylMatrix::PutVal Size" << Size(col);
		cout.flush();
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"Index out of range");
	} else if(index >= Size(col)) return 1;
	fElem[col][index] = value;
	//  delete[]newVet;
	this->fDecomposed = 0;
	return( 1 );
}

template<class TVar>
void TPZSkylMatrix<TVar>::MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
							const TVar alpha,const TVar beta ,const int opt,const int stride ) const {
	// Computes z = beta * y + alpha * opt(this)*x
	//          z and x cannot overlap in memory
	if ((!opt && this->Cols()*stride != x.Rows()) || this->Rows()*stride != x.Rows())
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__," <matrixs with incompatible dimensions>" );
	if(z.Rows() != x.Rows() || z.Cols() != x.Cols()) z.Redim(x.Rows(),x.Cols());
	if(x.Cols() != y.Cols() || x.Cols() != z.Cols() || x.Rows() != y.Rows() || x.Rows() != z.Rows()) {
		cout << "x.Cols = " << x.Cols() << " y.Cols()"<< y.Cols() << " z.Cols() " << z.Cols() << " x.Rows() " << x.Rows() << " y.Rows() "<< y.Rows() << " z.Rows() "<< z.Rows() << endl;
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__," incompatible dimensions\n");
	}
	this->PrepareZ(y,z,beta,opt,stride);
	int rows = this->Rows();
	int xcols = x.Cols();
	int ic, r;
	for (ic = 0; ic < xcols; ic++) {
		for( r = 0 ; r < rows ; r++ ) {
			int offset = Size(r);
			TVar val = 0.;
			const TVar *p = &x.g((r-offset+1)*stride,ic);
			TVar *diag = fElem[r] + offset-1;
			TVar *diaglast = fElem[r];
			while( diag > diaglast ) {
				val += *diag-- * *p;
				p += stride;
			}
			if( diag == diaglast ) val += *diag * *p;
			z(r*stride,ic) += val*alpha;
			TVar *zp = &z((r-offset+1)*stride,ic);
			val = x.g(r*stride,ic);
			diag = fElem[r] + offset-1;
			while( diag > diaglast ) {
				*zp += alpha * *diag-- * val;
				zp += stride;
			}
		}
	}
}

template<class TVar>
void TPZSkylMatrix<TVar>::SolveSOR(int & numiterations,const TPZFMatrix<TVar> &F,
							 TPZFMatrix<TVar> &result, TPZFMatrix<TVar> *residual, TPZFMatrix<TVar> &scratch,const TVar overrelax,
							 TVar &tol,const int FromCurrent,const int direction)  {
	
	if(residual == &F) {
		cout << "TPZMatrix::SolveSOR called with residual and F equal, no solution\n";
		return;
	}
	TVar res = 2*tol+1.;;
	if(residual) res = Norm(*residual);
	if(!FromCurrent) {
		result.Zero();
	}
	int r = this->Dim();
	int c = F.Cols();
	int i,ifirst = 0, ilast = r, iinc = 1;
	if(direction == -1) {
		ifirst = r-1;
		ilast = 0;
		iinc = -1;
	}
	int it;
	for(it=0; it<numiterations && res > tol; it++) {
		res = 0.;
		scratch = F;
		for(int ic=0; ic<c; ic++) {
			if(direction == 1) {
				//
				// compute the upper triangular part first and put into the scractch vector
				//
				for(i=ifirst; i!=ilast; i+= iinc) {
					//TPZColuna *mydiag = &fDiag[i];
					int offset = Size(i);
					TVar val;
					TVar *diag;
					TVar *diaglast = fElem[i];
					TVar *scratchp = &scratch(i-offset+1,ic);
					val = result(i,ic);
					diag = fElem[i] + offset-1;
					int lastid = diag-diaglast;
					int id;
					for(id=0; id<=lastid; id++) *(scratchp+id) -= *(diag-id) * val;
					/* codeguard fix
					 while( diag >= diaglast ) *scratchp++ -= *diag-- * val;
					 */
				}
				//
				// perform the SOR operation
				//
				for(i=ifirst; i!=ilast; i+= iinc) {
					//TPZColuna *mydiag = &fDiag[i];
					int offset = Size(i);
					TVar val = scratch(i,ic);
					TVar *p = &result(i-offset+1,ic);
					TVar *diag = fElem[i] + offset-1;
					TVar *diaglast = fElem[i];
					while( diag > diaglast ) val -= *diag-- * *p++;
					res += val*val;
					result(i,ic) += val*overrelax/ *diag;
				}
			} else {
				//
				// the direction is upward
				//
				// put the lower triangular part of the multiplication into the scratch vector
				//
				for(i=ifirst; i!=ilast; i+= iinc) {
					//TPZColuna *mydiag = &fDiag[i];
					int offset = Size(i);
					TVar val = scratch(i,ic);
					TVar *p = &result(i-offset+1,ic);
					TVar *diag = fElem[i] + offset-1;
					TVar *diaglast = fElem[i];
					while( diag > diaglast ) val -= *diag-- * *p++;
					//					res += val*val;
					scratch(i,ic) = val;
				}
				//
				// perform the SOR operation
				//
				for(i=ifirst; i!=ilast; i+= iinc) {
					//TPZColuna *mydiag = &fDiag[i];
					int offset = Size(i);
					//	REAL val = scratch(i,ic);
					TVar *diag;
					TVar *diaglast = fElem[i];
					TVar *scratchp = &scratch(i-offset+1,ic);
					//val= result(i,ic);
					TVar val = scratch(i,ic);
					val -= *diaglast * result(i,ic);
					res += val*val;
					val = overrelax * val / *diaglast;
					result(i,ic) += val;
					val = result(i,ic);
					diag = fElem[i] + offset-1;
					while( diag > diaglast ) *scratchp++ -= *diag-- * val;
				}
			}
		}
		res = sqrt(res);
	}
	if(residual) {
		this->Residual(result,F,*residual);
	}
	numiterations = it;
	tol = res;
}

/**************/
/*** GetVal ***/

template<class TVar>
const TVar &
TPZSkylMatrix<TVar>::GetVal(const int r,const int c ) const
{
	// inicializando row e col para trabalhar com a triangular superior
	int row(r),col(c);
	if ( row > col )
		this->Swap( &row, &col );
	
	if(row >= this->Dim() || col >= this->Dim() || row < 0 || col<0) {
		cout << "TPZSkylMatrix::GetVal index out of range row = " << row
		<< " col = " << col << endl;
		return this->gZero;
	}
	// Indice do vetor coluna.
	int index   = col - row;
	//TPZColuna *pCol = &fDiag[col];
	
	if ( index < Size(col) )
		return( fElem[col][index] );
	else {
		if(this->gZero != 0.) {
			cout << "TPZSkylMatrix gZero = " << this->gZero << endl;
			DebugStop();
		}
		return(this->gZero );
	}
}



/******** Operacoes com matrizes SKY LINE  ********/

/******************/
/*** Operator = ***/

template<class TVar>
TPZSkylMatrix<TVar> &
TPZSkylMatrix<TVar>::operator=(const TPZSkylMatrix<TVar> &A )
{
	Clear();
	Copy( A );
	return( *this );
}



/******************/
/*** Operator + ***/

template<class TVar>
TPZSkylMatrix<TVar>
TPZSkylMatrix<TVar>::operator+(const TPZSkylMatrix<TVar> &A) const
{
	if ( this->Dim() != A.Dim() )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"<incompatible dimensions>" );
	
	TPZVec<int> skylinesize(this->Dim());
	ComputeMaxSkyline(*this,A,skylinesize);
	TPZSkylMatrix res( this->fRow, skylinesize );
	
	TVar *elemMax;
	TVar *elemMin;
	int  sizeMax;
	int  sizeMin;
	
	for ( int col = 0; col < this->Dim(); col++ )
    {
		// Define o tamanho e os elementos da maior e da menor
		//  coluna.
		if ( Size(col) > A.Size(col) )
		{
			sizeMax = Size(col);
			elemMax = fElem[col];
			sizeMin = A.Size(col);
			elemMin = A.fElem[col];
		}
		else
		{
			sizeMax = A.Size(col);
			elemMax = A.fElem[col];
			sizeMin = Size(col);
			elemMin = fElem[col];
		}
		
		// Inicializa coluna da matriz resultado.
		
		// Efetua a SOMA.
		TVar *dest = res.fElem[col];
		int i;
		for ( i = 0; i < sizeMin; i++ )
			*dest++ = (*elemMax++) + (*elemMin++);
		for ( ; i < sizeMax; i++ )
			*dest++ = *elemMax++;
    }
	
	return( res );
}



/******************/
/*** Operator - ***/

template<class TVar>
TPZSkylMatrix<TVar>
TPZSkylMatrix<TVar>::operator-(const TPZSkylMatrix<TVar> &A ) const
{
	if ( this->Dim() != A.Dim() )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "operator-( TPZSkylMatrix ) <incompatible dimensions>" );
	
	TPZVec<int> skylinesize(this->Dim());
	ComputeMaxSkyline(*this,A,skylinesize);
	TPZSkylMatrix<TVar> res( this->fRow, skylinesize );
	
	for ( int col = 0; col < this->fRow; col++ )
    {
		// Define o tamanho e os elementos das colunas das 2 matrizes.
		int  sizeThis  = Size(col);
		TVar *elemThis = fElem[col];
		int  sizeA     = A.Size(col);
		TVar *elemA    = A.fElem[col];
		
		// Inicializa coluna da matriz resultado.
		
		// Efetua a SUBTRACAO.
		TVar *dest = res.fElem[col];
		int i;
		for ( i = 0; (i < sizeThis) && (i < sizeA); i++ ) *dest++ = (*elemThis++) - (*elemA++);
		if ( i == sizeA ) for ( ; i < sizeThis; i++ ) *dest++ = *elemThis++;
		else for ( ; i < sizeA; i++ ) *dest++ = -(*elemA++);
    }
	
	return( res );
}




/*******************/
/*** Operator += ***/

template<class TVar>
TPZSkylMatrix<TVar> &
TPZSkylMatrix<TVar>::operator+=(const TPZSkylMatrix<TVar> &A )
{
	if ( this->Dim() != A.Dim() )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"operator+=( TPZSkylMatrix ) <incompatible dimensions>" );
	
	TPZSkylMatrix res((*this)+A);
	*this = res;
	return *this;
}



/*******************/
/*** Operator -= ***/

template<class TVar>
TPZSkylMatrix<TVar> &
TPZSkylMatrix<TVar>::operator-=(const TPZSkylMatrix<TVar> &A )
{
	if ( this->Dim() != A.Dim() )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"operator-=( TPZSkylMatrix ) <incompatible dimensions>" );
	
	TPZSkylMatrix res(*this-A);
	*this = res;
	return *this;
}






/******** Operacoes com valores NUMERICOS ********/
//
// As operacoes com valores numericos sao efetuadas apenas nos
// elementos alocados. Em especial, as operacoes A = 0.0 e A *= 0.0
// desalocam todos os elementos da matriz.
//



/*****************************/
/*** Operator * ( REAL ) ***/

template<class TVar>
TPZSkylMatrix<TVar>
TPZSkylMatrix<TVar>::operator*(const TVar value ) const
{
	TPZSkylMatrix res( *this );
	
	for ( int col = 0; col < this->Dim(); col++ )
    {
		// Aloca nova coluna para o resultado.
		int colsize = Size(col);
		// Efetua a SOMA.
		TVar *elemRes  = res.fElem[col];
		for ( int i = 0; i < colsize; i++ )
			*elemRes++ *= value;
    }
	
	return( res );
}





/******************************/
/*** Operator *= ( REAL ) ***/

template<class TVar>
TPZSkylMatrix<TVar> &
TPZSkylMatrix<TVar>::operator*=(const TVar value )
{
	if ( IsZero( value ) )
    {
		Zero();
		return( *this );
    }
	
	int col, colmax = this->Dim();
	for (col=0; col<colmax; col++ )
    {
		// Efetua a MULTIPLICACAO.
		TVar *elem = fElem[col];
		TVar *end  = fElem[col+1];
		while ( elem < end ) *elem++ *= value;
    }
	
	this->fDecomposed = 0;
	return( *this );
}



/**************/
/*** Resize ***/
//
// Muda as dimensoes da matriz, mas matem seus valores antigos. Novas
// posicoes sao criadas com ZEROS.
//
template<class TVar>
int TPZSkylMatrix<TVar>::Resize( int newDim ,int ) {
	if ( newDim == this->Dim() )
		return( 1 );
	
	fElem.Resize(newDim+1);
	// Cria nova matrix.
	
	// Copia os elementos para a nova matriz.
	int min = MIN( newDim, this->Dim() );
	int i;
	for ( i = min+1; i <= newDim; i++ )
		fElem[i] = fElem[i-1];
	
	// Zera as posicoes que sobrarem (se sobrarem)
	fStorage.Resize(fElem[newDim]-fElem[0]);
	this->fRow = this->fCol = newDim;
	this->fDecomposed = 0;
	return( 1 );
}



/*************/
/*** Redim ***/
//
// Muda as dimensoes da matriz e ZERA seus elementos.
//
template<class TVar>
int
TPZSkylMatrix<TVar>::Redim( int newDim , int)
{
	if ( newDim == this->Dim() )
    {
		Zero();
		return( 1 );
    }
	
	Clear();
	fElem.Resize(newDim);
	fElem.Fill(0);
	this->fRow = this->fCol = newDim;
	this->fDecomposed = 0;
	return( 1 );
}



/**************************/
/*** Decompose Cholesky ***/
template<class TVar>
int
TPZSkylMatrix<TVar>::Decompose_Cholesky(std::list<int> &singular)
{
	if(this->fDecomposed == ECholesky) return 1;
	if (  this->fDecomposed )  TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Decompose_Cholesky <Matrix already Decomposed>" );
	
	singular.clear();
	TVar pivot;
	int dimension = this->Dim();
	/*  if(Dim() > 100) {
	 cout << "\nTPZSkylMatrix Cholesky decomposition Dim = " << Dim() << endl;
	 cout.flush();
	 }*/
	for ( int k = 0; k < dimension; k++ )
	{
		/*    if(!(k%100) && Dim() > 100) {
		 cout <<  k << ' ';
		 cout.flush();
		 }
		 if(!(k%1000)) cout << endl;*/
		if ( Size(k) == 0 )	return( 0 );
		
		// Faz sum = SOMA( A(k,p) * A(k,p) ), p = 1, ..., k-1.
		//
		TVar sum = 0.0;
		TVar *elem_k = fElem[k]+1;
		TVar *end_k  = fElem[k]+Size(k);
		for ( ; elem_k < end_k; elem_k++ ) sum += (*elem_k) * (*elem_k);
		
		// Faz A(k,k) = sqrt( A(k,k) - sum ).
		//
		pivot = fElem[k][0] - sum;
		if ( pivot < 1.e-10 ) {
			singular.push_back(k);
			pivot = 1.;
		}
		// A matriz nao e' definida positiva.
		
		pivot = fElem[k][0] = sqrt( pivot );
		
		// Loop para i = k+1 ... Dim().
		//
		int i=k+1;
		for ( int j = 2; i<dimension; j++,i++ ) {
			// Se tiverem elementos na linha 'i' cuja coluna e'
			//  menor do que 'K'...
			if ( Size(i) > j ) {
				// Faz sum = SOMA( A(i,p) * A(k,p) ), p = 1,..,k-1.
				sum = 0.0;
				TVar *elem_i = &fElem[i][j];
				TVar *end_i  = fElem[i+1];
				elem_k = &(fElem[k][1]);
				while ( (elem_i < end_i) && (elem_k < end_k) ) sum += (*elem_i++) * (*elem_k++);
				
				// Faz A(i,k) = (A(i,k) - sum) / A(k,k)
				fElem[i][j-1] = (fElem[i][j-1] -sum) / pivot;
			} else if ( Size(i) == j ) fElem[i][j-1] /= pivot;
			
			// Se nao tiverem estes elementos, sum = 0.0.
			
			// Se nao existir nem o elemento A(i,k), nao faz nada.
		}
	}
	
	if(this->Rows() && (GetVal(this->Rows()-1,this->Rows()-1)) < 1.e-15)
	{
		singular.push_back(this->Rows()-1);
		PutVal(this->Rows()-1,this->Rows()-1,1.);
	}
	this->fDecomposed  = ECholesky;
	this->fDefPositive = 1;
	return( 1 );
}

template<class TVar>
int
TPZSkylMatrix<TVar>::Decompose_Cholesky()
{
	if(this->fDecomposed == ECholesky) return 1;
	if (this->fDecomposed )  TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Decompose_Cholesky <Matrix already Decomposed>" );
	
	TVar pivot;
	int dimension = this->Dim();
	/*  if(Dim() > 100) {
	 cout << "\nTPZSkylMatrix Cholesky decomposition Dim = " << Dim() << endl;
	 cout.flush();
	 }*/
	for ( int k = 0; k < dimension; k++ )
    {
		/*      if(!(k%100) && Dim() > 100) {
		 cout <<  k << ' ';
		 cout.flush();
		 }
		 if(!(k%1000)) cout << endl;*/
		if ( Size(k) == 0 )	return( 0 );
		
		// Faz sum = SOMA( A(k,p) * A(k,p) ), p = 1, ..., k-1.
		//
		TVar sum = 0.0;
		TVar *elem_k = fElem[k]+1;
		TVar *end_k  = fElem[k]+Size(k);
		for ( ; elem_k < end_k; elem_k++ ) sum += (*elem_k) * (*elem_k);
		
		// Faz A(k,k) = sqrt( A(k,k) - sum ).
		//
		pivot = fElem[k][0] - sum;
		if ( pivot < 1.e-25 ) {
			cout << "TPZSkylMatrix::DecomposeCholesky a matrix nao e positiva definida" << pivot << endl;
			return( 0 );
		}
		// A matriz nao e' definida positiva.
		
		pivot = fElem[k][0] = sqrt( pivot );
		
		// Loop para i = k+1 ... Dim().
		//
		int i=k+1;
		for ( int j = 2; i<dimension; j++,i++ ) {
			// Se tiverem elementos na linha 'i' cuja coluna e'
			//  menor do que 'K'...
			if ( Size(i) > j ) {
				// Faz sum = SOMA( A(i,p) * A(k,p) ), p = 1,..,k-1.
				sum = 0.0;
				TVar *elem_i = &fElem[i][j];
				TVar *end_i  = fElem[i+1];
				elem_k = &(fElem[k][1]);
				while ( (elem_i < end_i) && (elem_k < end_k) ) sum += (*elem_i++) * (*elem_k++);
				
				// Faz A(i,k) = (A(i,k) - sum) / A(k,k)
				fElem[i][j-1] = (fElem[i][j-1] -sum) / pivot;
			} else if ( Size(i) == j ) fElem[i][j-1] /= pivot;
			
			// Se nao tiverem estes elementos, sum = 0.0.
			
			// Se nao existir nem o elemento A(i,k), nao faz nada.
		}
    }
	
	this->fDecomposed  = ECholesky;
	this->fDefPositive = 1;
	return( 1 );
}



/**********************/
/*** Decompose LDLt ***/
template<class TVar>
int
TPZSkylMatrix<TVar>::Decompose_LDLt(std::list<int> &singular)
{
	if( this->fDecomposed == ELDLt) return 1;
	if ( this->fDecomposed )
    {
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Decompose_LDLt <Matrix already Decomposed with different decomposition>" );
    }
	singular.clear();
	
	// Third try
	TVar *elj,*ell;
	int j,l,minj,minl,minrow,dimension = this->Dim();
	TVar sum;
	j = 1;
	while(j < dimension) {
		/*    if(!(j%100) && Dim() > 100) {
		 cout <<  j << ' ';
		 cout.flush();
		 }
		 if(!(j%1000)) cout << endl;*/
		minj = j-Size(j)+1;
		l = minj;
		while(l <= j) {
			minl = l-Size(l)+1;
			minrow = (minj<minl)? minl:minj;
			int k = minrow;
			//			DiagkPtr = fDiag+minrow;
			elj = fElem[j]+j-minrow;
			ell = fElem[l]+l-minrow;
			sum = 0.;
			while(k < l) {
				sum += *elj-- * *ell-- * *(fElem[k++]);
			}
			*elj -= sum;
			if(ell != elj) *elj /= *ell;
			else if(IsZero(*elj)) {
				singular.push_back(l);
				*elj = 1.;
			}
			l++;
		}
		j++;
	}
	
	if(this->Rows() && IsZero(GetVal(this->Rows()-1,this->Rows()-1)))
	{
		singular.push_back(this->Rows()-1);
		PutVal(this->Rows()-1,this->Rows()-1,1.);
	}
	this->fDecomposed  = ELDLt;
	this->fDefPositive = 0;
	//if(Dim() > 100) cout << endl;
	return( 1 );
}

template<class TVar>
int
TPZSkylMatrix<TVar>::Decompose_LDLt()
{
	
	if( this->fDecomposed == ELDLt) return 1;
	if (  this->fDecomposed )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Decompose_LDLt <Matrix already Decomposed with different decomposition>" );
	
	// Third try
	TVar *elj,*ell;
	int j,l,minj,minl,minrow,dimension = this->Dim();
	TPZVec<REAL> diag(dimension);
	for(j=0; j<dimension; j++)
	{
		diag[j] = *fElem[j];
	}
	TVar sum;
	j = 1;
	while(j < dimension) {
		/*    if(!(j%100) && Dim() > 100) {
		 cout <<  j << ' ';
		 cout.flush();
		 }
		 if(!(j%1000)) cout << endl;*/
		minj = j-Size(j)+1;
		l = minj;
		while(l <= j) {
			minl = l-Size(l)+1;
			minrow = (minj<minl)? minl:minj;
			int k = minrow;
			//			DiagkPtr = fDiag+minrow;
			elj = fElem[j]+j-minrow;
			ell = fElem[l]+l-minrow;
			TVar *diagptr = &diag[k];
			sum = 0.;
			while(k < l) {
				//		  sum += *elj-- * *ell-- * *(fElem[k++]);
				sum += *elj-- * *ell-- * *diagptr++;
				k++;
			}
			*elj -= sum;
			if(ell != elj) *elj /= *ell;
			else if(IsZero(*elj)) {
#ifdef LOG4CXX
				std::stringstream sout;
				sout << "col = " << j << " diagonal " << *elj;
				LOGPZ_DEBUG(logger,sout.str())
#endif
				
				*diagptr = *elj;
				cout << "TPZSkylMatrix pivot = " << *elj << endl;
				cout << "TPZSkylMatrix::DecomposeLDLt zero pivot\n";
				cout << "j = " << j << " l = " << l << endl;
			}
			else
			{
				*diagptr = *elj;
			}
			l++;
		}
		j++;
	}
	this->fDecomposed  = ELDLt;
	this->fDefPositive = 0;
	//if(Dim() > 100) cout << endl;
	return( 1 );
}



/*********************/
/*** Subst Forward ***/
//
//  Faz Ax = b, onde A e' triangular inferior.
//
template<class TVar>
int
TPZSkylMatrix<TVar>::Subst_Forward( TPZFMatrix<TVar> *B ) const
{
	if ( (B->Rows() != this->Dim()) || this->fDecomposed != ECholesky)
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"TPZSkylMatrix::Subst_Forward not decomposed with cholesky");
	
	//	std::cout << "SubstForward this " << (void *) this << " neq " << Dim() << " normb " << Norm(*B) << std::endl;
	int dimension=this->Dim();
    for ( int j = 0; j < B->Cols(); j++ )
	{
		int k=0;
		while (k<dimension && (*B)(k,j) == 0) {
			k++;
		}
		//		std::cout << "kstart " << k << std::endl;
		for (; k < dimension; k++ )
		{
			// Faz sum = SOMA( A[k,i] * B[i,j] ), para i = 1, ..., k-1.
			//
			TVar sum = 0.0;
			TVar *elem_ki = fElem[k]+1;
			TVar *end_ki  = fElem[k+1];
			TVar* BPtr = &(*B)(k,j);   //(k-1,j)
			//	for ( int i = k-1; elem_ki < end_ki; i-- )
			//	  sum += (*elem_ki++) * B->GetVal( i, j );
			
			while(elem_ki < end_ki) sum += (*elem_ki++) * (*--BPtr);//(*BPtr--)
			// Faz B[k,j] = (B[k,j] - sum) / A[k,k].
			//
			//	B->PutVal( k, j, (B->GetVal(k, j) - sum) / row_k->pElem[0] );
			BPtr = &(*B)(k,j);
			*BPtr-= sum;
			*BPtr /= fElem[k][0];
		}
	}
	
	return( 1 );
}

/*********************/
/*** Subst Backward ***/
//
//  Faz Ax = b, onde A e' triangular superior.
//
template<class TVar>
int
TPZSkylMatrix<TVar>::Subst_Backward( TPZFMatrix<TVar> *B ) const
{
	//	return TSimMatrix::Subst_Backward(B);
	
	if ( (B->Rows() != this->Dim()) || this->fDecomposed != ECholesky)
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"TPZSkylMatrix::Subst_Backward not decomposed with cholesky");
	
	//	std::cout << "SubstBackward this " << (void *) this << " neq " << Dim() << " ncols " << B->Cols() << std::endl;
	
	int Dimension = this->Dim();
	if(!Dimension) return 1;	// nothing to do
	int j;
    for ( j = 0; j < B->Cols(); j++ )
	{
		int k = Dimension-1;
		while (k>0 && (*B)(k,j) == 0.) {
			k--;
		}
		//		std::cout << "kstart " << k << std::endl;
		
		for (;k > 0; k-- )
		{
			// Faz sum = SOMA( A[k,i] * B[i,j] ), para i = 1, ..., k-1.
			//
			TVar val;
			TVar *elem_ki = fElem[k];
			TVar *end_ki  = fElem[k+1];
			TVar *BPtr = &(*B)(k,j);
			*BPtr /= *elem_ki++;
			val = *BPtr;
			//	BPtr;
			// substract the column of the skyline matrix from the vector.
			while(elem_ki < end_ki) *--BPtr -= (*elem_ki++) * val;
		}
	}
	for( j = 0; j< B->Cols(); j++) (*B)(0,j) /= fElem[0][0];
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
TPZSkylMatrix<TVar>::Subst_LForward( TPZFMatrix<TVar> *B ) const {
	if ( (B->Rows() != this->Dim()) || (this->fDecomposed != ELDLt && this->fDecomposed != ELU) )
		return( 0 );
	
	int dimension =this->Dim();
	for ( int k = 0; k < dimension; k++ ) {
		for ( int j = 0; j < B->Cols(); j++ ) {
			// Faz sum = SOMA( A[k,i] * B[i,j] ), para i = 1, ..., k-1.
			//
			TVar sum = 0.0;
			TVar *elem_ki = fElem[k]+1;
			TVar *end_ki  = fElem[k+1];
			TVar *BPtr = &(*B)(k,j);
			while(elem_ki < end_ki) sum += (*elem_ki++) * (*--BPtr);
			
			// Faz B[k,j] = (B[k,j] - sum) / A[k,k].
			//
			//	B->PutVal( k, j, (B->GetVal(k, j) - sum) / row_k->pElem[0] );
			BPtr = &(*B)(k,j);
			*BPtr-= sum;
		}
	}
	return( 1 );
}




/******************/
/*** Subst Diag ***/
//
//  Faz Ax = b, sendo que A e' assumida ser uma matriz diagonal.
//
template<class TVar>
int
TPZSkylMatrix<TVar>::Subst_Diag( TPZFMatrix<TVar> *B ) const {
	if ( (B->Rows() != this->Dim()) || this->fDecomposed != ELDLt) return( 0 );
	int dimension = this->Dim();
	for ( int j = 0; j < B->Cols(); j++ ) {
		TVar *BPtr = &(*B)(0,j);
		int k=0;
		while(k < dimension) *BPtr++ /= *(fElem[k++]);
	}
	return( 1 );
}

/*********************/
/*** Subst Backward ***/
//
//  Faz Ax = b, onde A e' triangular superior.
//
template<class TVar>
int
TPZSkylMatrix<TVar>::Subst_LBackward( TPZFMatrix<TVar> *B ) const
{
	//	return TSimMatrix::Subst_Backward(B);
	
	if ( (B->Rows() != this->Dim()) || !this->fDecomposed || this->fDecomposed == ECholesky)
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"TPZSkylMatrix::Subst_LBackward not decomposed properly");
	
	int Dimension = this->Dim();
	for ( int k = Dimension-1; k > 0; k-- ) {
		for ( int j = 0; j < B->Cols(); j++ ) {
			// Faz sum = SOMA( A[k,i] * B[i,j] ), para i = 1, ..., k-1.
			//
			//		REAL val = 0.0;
			TVar *elem_ki = fElem[k]+1;
			TVar *end_ki  = fElem[k+1];
			TVar *BPtr = &(*B)(k,j);
			// substract the column of the skyline matrix from the vector.
			
			TVar val = *BPtr;
			while(elem_ki < end_ki) *--BPtr -= (*elem_ki++) * val;
		}
	}
	return( 1 );
}





/**************************** PRIVATE ****************************/


/***************/
/*** PutZero ***/

template<class TVar>
int
TPZSkylMatrix<TVar>::Zero()
{
	
	fStorage.Fill(0.);
	this->fDecomposed = 0;
	this->fDefPositive = 0;
	return( 1 );
}



/*************/
/*** Error ***/
/*int
 TPZSkylMatrix::Error(const char *msg1,const char* msg2 )
 {
 ostringstream out;
 out << "TPZSkylMatrix::" << msg1 << msg2 << ".\n";
 //pzerror.Show();
 LOGPZ_ERROR (logger, out.str().c_str());
 DebugStop();
 return 0;
 }
 */


/*************/
/*** CLear ***/

template<class TVar>
int
TPZSkylMatrix<TVar>::Clear()
{
	this->fStorage.Resize(0);
	fStorage.Shrink();
	this->fElem.Resize(0);
	this->fRow = this->fCol = 0;
	this->fDecomposed = 0;
	return( 1 );
}



/************/
/*** Copy ***/

template<class TVar>
void
TPZSkylMatrix<TVar>::Copy(const TPZSkylMatrix<TVar> &A )
{
	int dimension = A.Dim();
	this->fRow = this->fCol = dimension;
	fElem.Resize(A.fElem.NElements());
	fStorage = A.fStorage;
	int i;
	TVar *firstp = 0;
	if(fStorage.NElements()) firstp = &fStorage[0];
	for(i=0; i<=dimension; i++)
		fElem[i]=firstp+(A.fElem[i]-A.fElem[0]);
	this->fDecomposed  = A.fDecomposed;
	this->fDefPositive = A.fDefPositive;
	
}

template<class TVar>
TPZSkylMatrix<TVar>
TPZSkylMatrix<TVar>::operator-() const { return operator*(-1.0); }


#ifdef OOPARLIB

int TPZSkylMatrix<TVar>::Unpack (TReceiveStorage *buf ){
	TSimMatrix::Unpack(buf);
	int rows;
	buf->UpkInt(&rows);
	Redim(rows,rows);
	int sz;
	for(int i=0;i<rows;i++) {
		sz=fDiag[i].size;
		buf->UpkInt(&sz);
		for (int j=0;j<sz;j++) {
			buf->UpkDouble(fDiag[i].pElem);
		}
	}
	return 1;
}

int TPZSkylMatrix::Pack( TSendStorage *buf ){
	TSimMatrix::Pack(buf);
	int rows = Rows();
	int sz;
	for (int i=0;i<rows;i++) {
		buf->PkInt(&sz);
		fDiag[i].size=sz;
		fDiag[i].vetSize=sz;
		fDiag[i].pElem=new REAL[sz];
		for (int j=0;j<sz;j++) {
			buf->PkDouble(fDiag[i].pElem);
		}
	}
	return 1;
}

TSaveable *TPZSkylMatrix::Restore(TReceiveStorage *buf) {
	TPZSkylMatrix *m = new TPZSkylMatrix(0);
	m->Unpack(buf);
	return m;
}


int TPZSkylMatrix::DerivedFrom(long Classid){
	if(Classid == GetClassID()) return 1;
	return TSimMatrix::DerivedFrom(Classid);
}

int TPZSkylMatrix::DerivedFrom(char *classname){
	if(!strcmp(ClassName(),classname)) return 1;
	return TSimMatrix::DerivedFrom(classname);
}

#endif

template<class TVar>
void TPZSkylMatrix<TVar>::DecomposeColumn(int col, int prevcol){
	TVar *ptrprev;     //Pointer to prev column
	TVar *ptrcol;      //Pointer to col column
	int skprev, skcol; //prev and col Skyline height respectively
	int minline;
	
	skprev = SkyHeight(prevcol);
	skcol = SkyHeight(col);
	
	ptrprev = Diag(prevcol);
	ptrcol = Diag(col);
	
	if((prevcol-skprev) > (col-skcol)){
		minline = prevcol - skprev;
	}else
    {
		minline = col - skcol;
    }
	if(minline > prevcol) {
		cout << "error condition\n";
		cout.flush();
		return;
	}
	TVar *run1 = ptrprev + (prevcol-minline);
	TVar *run2 = ptrcol + (col-minline);
	TVar sum = 0;
	/*
	 while(run1-ptrprev > templatedepth) {
	 run1-=templatedepth-1;
	 run2-=templatedepth-1;
	 sum += TemplateSum<templatedepth>(run1--,run2--);
	 }
	 */
	
	while(run1 != ptrprev) sum += (*run1--)*(*run2--);
	*run2-=sum;
	if(run1 != run2){
		*run2 /= *run1;
	}else{
		*run2=sqrt(*run2);
	}
	
}

template<class TVar>
void TPZSkylMatrix<TVar>::DecomposeColumn(int col, int prevcol,std::list<int> &singular){
	TVar *ptrprev;     //Pointer to prev column
	TVar *ptrcol;      //Pointer to col column
	int skprev, skcol; //prev and col Skyline height respectively
	int minline;
	
	skprev = SkyHeight(prevcol);
	skcol = SkyHeight(col);
	
	ptrprev = Diag(prevcol);
	ptrcol = Diag(col);
	
	if((prevcol-skprev) > (col-skcol)){
		minline = prevcol - skprev;
	}else
    {
		minline = col - skcol;
    }
	if(minline > prevcol) {
		cout << "error condition\n";
		cout.flush();
		return;
	}
	TVar *run1 = ptrprev + (prevcol-minline);
	TVar *run2 = ptrcol + (col-minline);
	TVar sum = 0;
	/*
	 while(run1-ptrprev > templatedepth) {
	 run1-=templatedepth-1;
	 run2-=templatedepth-1;
	 sum += TemplateSum<templatedepth>(run1--,run2--);
	 }
	 */
	
	while(run1 != ptrprev) sum += (*run1--)*(*run2--);
	*run2-=sum;
	if(run1 != run2){
		*run2 /= *run1;
	}else{
		TVar pivot = *run2;
		if ( pivot < 1.e-10 ) {
#ifdef LOG4CXX
			std::stringstream sout;
			sout << "equation " << col << " is singular pivot " << pivot;
			LOGPZ_WARN(logger,sout.str())
#endif
			singular.push_back(col);
			pivot = 1.;
		}
		
		*run2=sqrt(pivot);
	}
	
}

template<class TVar>
void TPZSkylMatrix<TVar>::DecomposeColumn2(int col, int prevcol){
	
	//cout << "lcol " << lcol << " with " << lprevcol << endl;
	//cout.flush();
	//Cholesky Decomposition
	TVar *ptrprev;     //Pointer to prev column
	TVar *ptrcol;      //Pointer to col column
	int skprev, skcol; //prev and col Skyline height respectively
	int minline;
	
	skprev = SkyHeight(prevcol);
	skcol = SkyHeight(col);
	
	ptrprev = Diag(prevcol);
	ptrcol = Diag(col);
	
	if((prevcol-skprev) > (col-skcol)){
		minline = prevcol - skprev;
	}else
    {
		minline = col - skcol;
    }
	if(minline > prevcol) {
		cout << "error condition\n";
		cout.flush();
		return;
	}
	TVar *run1 = ptrprev + 1;
	TVar *run2 = ptrcol + (col-prevcol)+1;
	TVar *lastptr = ptrprev + prevcol-minline+1;
	TVar sum = 0;
	TVar *modify = ptrcol+(col-prevcol);
#ifndef BLAS
	/*
	 while(lastptr-run1 > templatedepth) {
	 sum += TemplateSum<templatedepth>(run1,run2);
	 run1+=templatedepth;
	 run2+=templatedepth;
	 }
	 */
	
	while(run1 != lastptr) sum += (*run1++)*(*run2++);
	
	//cout << "col " << col << " prevcol " << prevcol << " sum " << sum << " modify " << *modify << " ptrprev " << *ptrprev;
#else
	int n=lastptr-run1;
	sum = cblas_ddot(n,run1,1,run2,1);
#endif
	*modify-=sum;
	if(col != prevcol){
		*modify /= *ptrprev;
	}else{
		if ( *modify < 1.e-25 ) {
			cout << "TPZSkylMatrix::DecomposeCholesky a matrix nao e positiva definida" << *modify << endl;
			*modify = 1.e-10;
			//      return;
		}
		
		*modify=sqrt(*modify);
	}
	//cout << " modified " << *modify << endl;
	//cout.flush();
	
}

/*
 template<class TVar>
 void TPZSkylMatrix<TVar>::TestSpeed(int col, int prevcol){
 
 TVar *ptrprev;     //Pointer to prev column
 TVar *ptrcol;      //Pointer to col column
 int skprev, skcol; //prev and col Skyline height respectively
 int minline;
 
 skprev = SkyHeight(prevcol);
 skcol = SkyHeight(col);
 
 ptrprev = Diag(prevcol);
 ptrcol = Diag(col);
 
 if((prevcol-skprev) > (col-skcol)){
 minline = prevcol - skprev;
 }else
 {
 minline = col - skcol;
 }
 if(minline > prevcol) {
 cout << "error condition\n";
 cout.flush();
 return;
 }
 int i, j;
 struct timeb meas1,meas2;
 ftime(&meas1);
 for(i=0;i<1000;i++) {
 for(j=0;j<1000;j++) {
 TVar *run1 = ptrprev + (prevcol-minline);
 TVar *run2 = ptrcol + (col-minline);
 TVar
 sum = 0;
 
 
 while(run1-ptrprev > templatedepth) {
 run1-=templatedepth-1;
 run2-=templatedepth-1;
 sum += TemplateSum<templatedepth>(run1--,run2--);
 }
 
 
 
 while(run1 != ptrprev) sum += (*run1--)*(*run2--);
 *run2-=sum;
 if(*run1 != *run2){
 *run2 /= *run1;
 }else{
 *run2=sqrt(*run2);
 }
 }
 }
 ftime(&meas2);
 int dif = (meas2.time-meas1.time)*1000 + meas2.millitm-meas1.millitm;
 cout << "Time elapsed in milliseconds " << dif << endl;
 
 }
 */

template class TPZSkylMatrix<REAL>;
