/**
 * @file
 * @brief Contains the implementation of the TPZBlockDiagonal methods.
 */
//
// Author: PHILIPPE DEVLOO
//
// File:   pzblockdiag.c
//
// Class:  TPZBlockDiagonal
//
// Obs.:   Implements block diagonal matrices
//
// Versao: 01 - 2001
//
#include "pzfmatrix.h"
#include "pzblockdiag.h"

//#include "pzerror.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "pzlog.h"
#include <sstream>
#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.StrMatrix"));
#endif

using namespace std;

/** @brief Initializing variable zero as real */
static REAL zero = 0.;

template<class TVar>

void TPZBlockDiagonal<TVar>::AddBlock(int i, TPZFMatrix<TVar> &block){
	// ::cout << "Iniciando inserc� de bloco na posic�\t" ;
	int firstpos = fBlockPos[i];
	//  ::cout << firstpos << "\n";
	int bsize = fBlockSize[i];
	//  ::cout << "Dimens� do bloco a ser inserido\t" << bsize << "\n";
	
	int r,c;
	for(r=0; r<bsize; r++) {
		for(c=0; c<bsize; c++) {
			//  ::cout << " inserindo elemento: \t local[" << r <<"," << c
			//     << "]\t global [" << firstpos+r+bsize*c << "]= \t"
			//     << block(r,c) << "\t totalizando \t" <<fStorage[firstpos+r+bsize*c] + block(r,c)
			//     << "\n";
			fStorage[firstpos+r+bsize*c] += block(r,c);
		}
	}
}

template<class TVar>
void TPZBlockDiagonal<TVar>::SetBlock(int i, TPZFMatrix<TVar> &block){
	// ::cout << "Iniciando inserc� de bloco na posic�\t" ;
	int firstpos = fBlockPos[i];
	//  ::cout << firstpos << "\n";
	int bsize = fBlockSize[i];
	//  ::cout << "Dimens� do bloco a ser inserido\t" << bsize << "\n";
	
	int r,c;
	for(r=0; r<bsize; r++) {
		for(c=0; c<bsize; c++) {
			//  ::cout << " inserindo elemento: \t local[" << r <<"," << c
			//     << "]\t global [" << firstpos+r+bsize*c << "]= \t"
			//     << block(r,c) << "\t totalizando \t" <<fStorage[firstpos+r+bsize*c] + block(r,c)
			//     << "\n";
			fStorage[firstpos+r+bsize*c] = block(r,c);
		}
	}
}

template<class TVar>
void TPZBlockDiagonal<TVar>::GetBlock(int i, TPZFMatrix<TVar> &block){
    int firstpos = fBlockPos[i];
    int bsize = fBlockSize[i];
    block.Redim(bsize,bsize);
    int r,c;
    for(r=0; r<bsize; r++) {
        for(c=0; c<bsize; c++) {
            block(r,c) = fStorage[firstpos+r+bsize*c];
        }
    }
}
template<class TVar>
void TPZBlockDiagonal<TVar>::Initialize(const TPZVec<int> &blocksize){
	int nblock = blocksize.NElements();
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "Number of blocks \t" << nblock;
		LOGPZ_DEBUG(logger,sout.str());
	}
#endif
	//  ::cout << "Nmero de blocos \t" << nblock << "\n";
	fBlockSize = blocksize;
	fBlockPos.Resize(nblock+1,0); 
	int b;
	int ndata = 0;
	int neq = 0;
	int bsize;
	for(b=0; b<nblock; b++) {
		bsize = blocksize[b];
		//    ::cout << "Dimens� do bloco [" << b << "]\t" <<bsize <<"\n";
		fBlockPos[b+1] = fBlockPos[b]+bsize*bsize;
		ndata += bsize*bsize;
		neq += bsize;
		//    ::cout << "Nmero de equac�s\t" << neq << "\n";
	}
#ifdef LOG4CXX
	if(ndata > 10000000)
	{
		std::stringstream sout;
		sout << "Calling fStorage.Resize(ndata,0.) with ndata = " << ndata;
		LOGPZ_DEBUG(logger,sout.str());
	}
#endif
	
	fStorage.Fill(0.,0);
	fStorage.Resize(ndata,0.);
	this->fDecomposed = 0;
	this->fRow = neq;
	this->fCol = neq;
	
	//fStorage.Print(std::cout);
}

template<class TVar>
void TPZBlockDiagonal<TVar>::BuildFromMatrix(TPZMatrix<TVar> &mat) {
	if(mat.Rows() != this->Rows()) {
		cout << "TPZBlockDiagonal::BuildFromMatrix wrong data structure\n";
		return;
	}
	int nblock = fBlockSize.NElements();
	int b,eq=0;
	for(b=0; b<nblock; b++) {
		int r,c,bsize = fBlockSize[b];
		int pos = fBlockPos[b];
		for(r=0; r<bsize; r++){
			for(c=0; c<bsize; c++) {
				fStorage[pos+r+c*bsize] = mat.GetVal(eq+r,eq+c);
			}
		}
		eq += bsize;
	}
}

/*******************/
/*** Constructor ***/

template<class TVar>
TPZBlockDiagonal<TVar>::TPZBlockDiagonal()
: TPZMatrix<TVar>(), fStorage(), fBlockPos(1,0),fBlockSize()
{
	
}

template<class TVar>
TPZBlockDiagonal<TVar>::TPZBlockDiagonal(const TPZVec<int> &blocksize)
: TPZMatrix<TVar>(), fStorage(), fBlockPos(1,0),fBlockSize()
{
	Initialize(blocksize);
}



/********************/
/*** Constructors ***/

template<class TVar>
TPZBlockDiagonal<TVar>::TPZBlockDiagonal(const TPZVec<int> &blocksizes, const TPZFMatrix<TVar> &glob)
: TPZMatrix<TVar>(), fBlockSize(blocksizes)
{
	int nblock = blocksizes.NElements();
	fBlockPos.Resize(nblock+1,0);
	int b;
	int ndata = 0;
	int neq = 0;
	int bsize;
	for(b=0; b<nblock; b++) {
		bsize = blocksizes[b];
		fBlockPos[b+1] = fBlockPos[b]+bsize*bsize;
		ndata += bsize*bsize;
		neq += bsize;
	}
	fStorage.Resize(ndata,0.);
	this->fRow = neq;
	this->fCol = neq;
	int pos;
	int eq = 0, r, c;
	for(b=0; b<nblock; b++) {
		bsize = fBlockSize[b];
		pos = fBlockPos[b];
		for(r=0; r<bsize; r++) {
			for(c=0; c<bsize; c++) {
				fStorage[pos+r+bsize*c]= glob.GetVal(eq+r,eq+c);
			}
		}
		eq += bsize;
	}
	
}



/*********************************/
/*** Constructor( TPZBlockDiagonal& ) ***/
template<class TVar>
TPZBlockDiagonal<TVar>::TPZBlockDiagonal (const TPZBlockDiagonal<TVar> & A)
: TPZMatrix<TVar>( A.Dim(), A.Dim() ), fStorage(A.fStorage),
fBlockPos(A.fBlockPos), fBlockSize(A.fBlockSize)
{
}



/******************/
/*** Destructor ***/

template<class TVar>
TPZBlockDiagonal<TVar>::~TPZBlockDiagonal ()
{
}



/***********/
/*** Put ***/
template<class TVar>
int
TPZBlockDiagonal<TVar>::Put(const int row,const int col,const TVar& value )
{
	//  cout << "TPZBlockDiagonal.Put should not be called\n";
	if ( (row >= Dim()) || (col >= Dim()) || row<0 || col<0)
    {
		cout << "TPZBlockDiagonal::Put: " << row << "," << col << "," << Dim();
		cout << "\n";
		return( 0 );
    }
	
	return( PutVal( row, col, value ) );
}

/***********/
/*** PutVal ***/
template<class TVar>
int
TPZBlockDiagonal<TVar>::PutVal(const int row,const int col,const TVar& value )
{
	//  cout << "TPZBlockDiagonal.PutVal should not be called\n";
	
	int b = 0;
	int nb = fBlockSize.NElements();
	if(nb==0) {
		cout << "TPZBlockDiagonal::PutVal called with parameters out of range\n";
		return -1;
	}
	int eq=0;
	int bsize = fBlockSize[b];
	while(eq+bsize <= row && b < nb) {
		eq+=bsize;
		b++;
		bsize = fBlockSize[b];
	}
	if(b==nb) {
    	cout << "TPZBlockDiagonal::PutVal wrong data structure\n";
		return -1;
	}
	if(col < eq || col >= eq+bsize) {
		if(value != 0.) {
    		cout << "TPZBlockDiagonal::PutVal, indices row col out of range\n";
			return -1;
		} else {
			return 0;
		}
	}
	fStorage[fBlockPos[b]+row-eq+bsize*(col-eq)] = value;
	return 0;
}



/***********/
/*** Get ***/

template<class TVar>
const TVar&
TPZBlockDiagonal<TVar>::Get(const int row,const int col ) const
{
	//  cout << "TPZBlockDiagonal::Get should never be called\n";
	if ( (row >= Dim()) || (col >= Dim()) )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "TPZBlockDiagonal::Get <indices out of band matrix range>" );
	
	return( GetVal( row, col ) );
}

template<class TVar>
TVar &
TPZBlockDiagonal<TVar>::operator()(const int row, const int col) {
	//  cout << "TPZBlockDiagonal.operator() should not be called\n";
	
	int b = 0;
	int nb = fBlockSize.NElements();
	if(nb==0) {
		cout << "TPZBlockDiagonal::operator() called with parameters out of range\n";
		zero = 0.;
		return zero;
	}
	int eq=0;
	int bsize = fBlockSize[b];
	while(eq+bsize <= row && b < nb) {
		eq+=bsize;
		b++;
		bsize = fBlockSize[b];
	}
	if(b==nb) {
		cout << "TPZBlockDiagonal::operator() wrong data structure\n";
		zero = 0.;
		return zero;
	}
	if(col < eq || col >= eq+bsize) {
		cout << "TPZBlockDiagonal::operator(), indices row col out of range\n";
		zero = 0.;
		return zero;
	}
	return fStorage[fBlockPos[b]+row-eq+bsize*(col-eq)];
}

/***********/
/*** GetVal ***/

template<class TVar>
const TVar &
TPZBlockDiagonal<TVar>::GetVal(const int row,const int col ) const
{
	// cout << "TPZBlockDiagonal.GetVal should not be called\n";
	
	int b = 0;
	int nb = fBlockSize.NElements();
	if(nb==0) {
		cout << "TPZBlockDiagonal::GetVal called with parameters out of range\n";
	}
	int eq=0;
	int bsize = fBlockSize[b];
	while(eq+bsize <= row && b < nb) {
		eq+=bsize;
		b++;
		bsize = fBlockSize[b];
	}
	if(b==nb) {
		cout << "TPZBlockDiagonal::GetVal wrong data structure\n";
	}
	if(col < eq || col >= eq+bsize) {
		//cout << "TPZBlockDiagonal::GetVal, indices row col out of range\n";
		zero = 0.;
		return zero;
	}
	return fStorage[fBlockPos[b]+row-eq+bsize*(col-eq)];
}




/******** Operacoes com MATRIZES GENERICAS ********/


/*******************/
/*** MultiplyAdd ***/
//
//  perform a multiply add operation to be used by iterative solvers
//

template<class TVar>
void TPZBlockDiagonal<TVar>::MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
							   const TVar alpha,const TVar beta ,const int opt,const int stride ) const {
	// Computes z = beta * y + alpha * opt(this)*x
	//          z and x cannot overlap in memory
	
	if ((!opt && this->Cols()*stride != x.Rows()) || this->Rows()*stride != x.Rows())
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "TPZBlockDiagonal::MultAdd <matrixs with incompatible dimensions>" );
	if(x.Cols() != y.Cols() || x.Cols() != z.Cols() || x.Rows() != y.Rows() || x.Rows() != z.Rows()) {
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"TPZBlockDiagonal::MultAdd incompatible dimensions\n");
	}
	
	this->PrepareZ(y,z,beta,opt,stride);
	//  int rows = Rows();
	int xcols = x.Cols();
	int nb= fBlockSize.NElements();
	int ic, b, bsize, eq=0, r, c;
	if(opt == 0) {
		for (ic = 0; ic < xcols; ic++) {
			eq=0;
			for(b=0; b<nb; b++) {
				bsize = fBlockSize[b];
				int pos = fBlockPos[b];
				for(r=0; r<bsize; r++) {
					for(c=0; c<bsize; c++) {
						z(eq+r,ic) += alpha*fStorage[pos+r+bsize*c]*x.GetVal((eq+c)*stride,ic);
					}
				}
				eq += bsize;
			}
		}
	} else {
		cout << "xcols \t" << xcols << "\n";
		for (ic = 0; ic < xcols; ic++) {
			eq=0;
			for(b=0; b<nb; b++) {
				bsize = fBlockSize[b];
				int pos = fBlockPos[b];
				for(r=0; r<bsize; r++) {
					for(c=0; c<bsize; c++) {
						z(eq+r,ic) += alpha*fStorage[pos+r+bsize*c]*x.GetVal((eq+c)*stride,ic);
						//   ::cout << "Z[" << (eq+r) <<"," << ic <<"] = " <<z(eq+r,ic) <<"\n";  
					}
				}
				eq+=bsize;
			}
		}
	}
}







/***************/
/**** Zero ****/
template<class TVar>
int
TPZBlockDiagonal<TVar>::Zero()
{
	
	fStorage.Fill(0.,0);
	this->fDecomposed = 0;
	
	return( 1 );
}



/********************/
/*** Transpose () ***/
template<class TVar>
void
TPZBlockDiagonal<TVar>::Transpose (TPZMatrix<TVar> *const T) const
{
	T->Resize( Dim(), Dim() );
	
	int b, bsize, eq = 0, pos;
	int nb = fBlockSize.NElements(), r, c;
	for ( b=0; b<nb; b++) {
		pos= fBlockPos[b];
		bsize = fBlockSize[b];
		for(r=0; r<bsize; r++) {
			for(c=0; c<bsize; c++) {
				T->PutVal(eq+r,eq+c,fStorage[pos+c+r*bsize]);
			}
		}
		eq += bsize;
	}
}


/*****************/
/*** Decompose_LU ***/
//fElem[ fBand * (2*row + 1) + col ]
template<class TVar>
int
TPZBlockDiagonal<TVar>::Decompose_LU(std::list<int> &singular)
{
	return Decompose_LU();
}

template<class TVar>
int
TPZBlockDiagonal<TVar>::Decompose_LU()
{
	
	LOGPZ_DEBUG(logger, "TPZBlockDiagonal::Decompose_LU");
	
	if (  this->fDecomposed && this->fDecomposed == ELU) {
		return ELU;
	} else if(this->fDecomposed) {
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"TPZBlockDiagonal::Decompose_LU is already decomposed with other scheme");
	}
	
	int b,nb,pos,bsize;
	nb = fBlockSize.NElements();
	for(b=0;b<nb; b++) {
		pos = fBlockPos[b];
		bsize = fBlockSize[b];
		if(!bsize) continue;
		
#ifdef DEBUG
		std::stringstream mess;
		mess << "TPZBlockDiagonal::Decompose_LU() - bsize = " << bsize << ", bsize*bsize = " << bsize*bsize;
		LOGPZ_DEBUG(logger,mess.str());
#endif
		
		TPZFMatrix<TVar> temp(bsize,bsize,&fStorage[pos],bsize*bsize);
		std::list<int> singular;
		temp.Decompose_LU(singular);
	}
	this->fDecomposed = ELU;
	return 1;
}

template<class TVar>
int
TPZBlockDiagonal<TVar>::Substitution( TPZFMatrix<TVar> *B) const
{
	if(this->fDecomposed != ELU) {
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"TPZBlockDiagonal::Decompose_LU is decomposed with other scheme");
	}
	
	int b,nb,pos,bsize,eq=0;
	nb = fBlockSize.NElements();
	int c, nc = B->Cols();
	for(c=0; c<nc; c++) {
		eq = 0;
		for(b=0;b<nb; b++) {
			pos = fBlockPos[b];
			bsize = fBlockSize[b];
			if(!bsize) continue;
			//      TPZFMatrix<>temp(bsize,bsize,&fStorage[pos],bsize*bsize);
			//      temp.SetIsDecomposed(ELU);
			TPZFMatrix<TVar> BTemp(bsize,1,&(B->operator()(eq,c)),bsize);
			TPZFMatrix<TVar>::Substitution(fStorage+pos,bsize,&BTemp);
			eq+= bsize;
		}
	}
	return 1;
}

/************************** Private **************************/

/*************/
/*** Error ***/

/*int
 TPZBlockDiagonal::Error(const char *msg1,const char *msg2 ) 
 {
 ostringstream out;
 out << "TPZBlockDiagonal::" << msg1 << msg2 << ".\n";
 LOGPZ_ERROR (logger, out.str().c_str());
 // pzerror.Show();
 //DebugStop();
 return 0;
 }*/



/*************/
/*** Clear ***/
template<class TVar>
int
TPZBlockDiagonal<TVar>::Clear()
{
	fStorage.Resize(0);
	fBlockPos.Resize(0);
	fBlockSize.Resize(0);
	this->fRow = 0;
	this->fCol = 0;
	this->fDecomposed = 0;
	return( 1 );
}

template<class TVar>
int TPZBlockDiagonal<TVar>::main() {
	
	cout << "Entering the main program\n";
	cout.flush();
	TPZFMatrix<TVar> ref(7,7,0.);
	int r,c;
	for(r=0; r<7; r++) {
		for(c=0; c<7; c++) {
			ref(r,c) = 5+r*c+3*r;
		}
		ref(r,r) += 1000;
	}
	TPZVec<int> blocksize(3);
	blocksize[0] = 2;
	blocksize[1] = 4;
	blocksize[2] = 1;
	TPZBlockDiagonal bd1(blocksize,ref);
	TPZBlockDiagonal bd2(bd1);
	ref.Print("original matrix");
	bd1.Solve_LU(&ref);
	bd1.Solve_LU(&ref);
	ref.Print("after inverting the diagonal");
	TPZFMatrix<TVar> ref2;
	bd2.Multiply(ref,ref2);
	bd2.Multiply(ref2,ref);
	ref.Print("restoring the original matrix");
	return 1;
	
}

template<class TVar>
void TPZBlockDiagonal<TVar>::Print(const char *msg, std::ostream &out, const MatrixOutputFormat format) const {
	
	if(format != EFormatted)
	{
		TPZMatrix<TVar>::Print(msg,out,format);
		return;
	}
	out << "TPZBlockDiagonal matrix ";
	if(msg) out << msg;
	out  << std::endl;
	
	int nblock = fBlockSize.NElements();
	out << "Number of blocks " << nblock << std::endl; 
	int b,bsize,pos;
	for(b=0; b<nblock; b++) {
		bsize = fBlockSize[b];
		out << "block number " << b << " size : " << bsize << std::endl;
		int r,c;
		pos = fBlockPos[b];
		for(c=0; c<bsize; c++) {
			for(r=0; r<bsize ; r++) {
				out << fStorage[pos+r+bsize*c] << ' ';
			}
			out << std::endl;
		}
	}
}

/**
 * Updates the values of the matrix based on the values of the matrix
 */
template<class TVar>
void TPZBlockDiagonal<TVar>::UpdateFrom(TPZAutoPointer<TPZMatrix<TVar> > mat)
{
	if(!mat) 
	{
		cout << "TPZBlockDiagonal::UpdateFrom" << " called with zero argument\n";
		return;
	}
	this->fDecomposed = ENoDecompose;
	int nblock = fBlockSize.NElements();
	int b,bsize,pos,firsteq = 0;
	for(b=0; b<nblock; b++) {
		bsize = fBlockSize[b];
		//    int r,c;
		pos = fBlockPos[b];
		TPZFMatrix<TVar> block(bsize,bsize,&fStorage[pos],bsize*bsize);
		mat->GetSub(firsteq,firsteq,bsize,bsize,block);
		firsteq += bsize;
	}
}

/** Fill the matrix with random values (non singular matrix) */
template<class TVar>
void TPZBlockDiagonal<TVar>::AutoFill() {

	int b, bsize, eq = 0, pos;
	int nb = fBlockSize.NElements(), r, c;
	for ( b=0; b<nb; b++) {
		pos= fBlockPos[b];
		bsize = fBlockSize[b];
		for(r=0; r<bsize; r++) {
            REAL sum = 0.;
			for(c=0; c<bsize; c++) {
                REAL val = ((REAL)rand())/RAND_MAX;
				fStorage[pos+c+r*bsize] = val;
                if(c!= r) sum += fabs(val);
			}
            if (fabs(fStorage[pos+r+r*bsize]) < sum) {
                fStorage[pos+r+r*bsize] = sum+1.;
            }
		}
		eq += bsize;
	}
}

template class TPZBlockDiagonal<REAL>;