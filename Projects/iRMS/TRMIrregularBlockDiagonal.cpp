//
//  TRMIrregularBlockDiagonal.cpp
//  PZ
//
//  Created by Omar on 5/26/16.
//
//

#include "TRMIrregularBlockDiagonal.h"


#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.StrMatrix"));
#endif

using namespace std;

template<class TVar>
void TRMIrregularBlockDiagonal<TVar>::AddBlock(long i, TPZFMatrix<TVar> &block){
    
    long firstpos = fBlockPos[i];
    long b_isize = fBlockSize[i].first;
    long b_jsize = fBlockSize[i].second;
    
    long r,c;
    for(r=0; r<b_isize; r++) {
        for(c=0; c<b_jsize; c++) {
            fStorage[firstpos+r+b_isize*c] += block(r,c);
        }
    }
}

template<class TVar>
void TRMIrregularBlockDiagonal<TVar>::SetBlock(long i, TPZFMatrix<TVar> &block){
    
    long firstpos = fBlockPos[i];
    long b_isize = fBlockSize[i].first;
    long b_jsize = fBlockSize[i].second;
    
    long r,c;
    for(r=0; r<b_isize; r++) {
        for(c=0; c<b_jsize; c++) {
            fStorage[firstpos+r+b_isize*c] = block(r,c);
        }
    }
}

template<class TVar>
void TRMIrregularBlockDiagonal<TVar>::GetBlock(long i, TPZFMatrix<TVar> &block){

    long firstpos = fBlockPos[i];
    long b_isize = fBlockSize[i].first;
    long b_jsize = fBlockSize[i].second;
    
    block.Redim(b_isize,b_jsize);
    
    long r,c;
    for(r=0; r<b_isize; r++) {
        for(c=0; c<b_jsize; c++) {
            block(r,c) = fStorage[firstpos+r+b_isize*c];
        }
    }
}
template<class TVar>
void TRMIrregularBlockDiagonal<TVar>::Initialize(TPZVec< std::pair<long, long> > &blocksize){
    long nblock = blocksize.NElements();
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << "Number of blocks \t" << nblock;
        LOGPZ_DEBUG(logger,sout.str());
    }
#endif
    fBlockSize = blocksize;
    fBlockPos.Resize(nblock+1,0);
    long b;
    long ndata = 0;
    long nr = 0, nc = 0;
    int b_isize, b_jsize;
    for(b=0; b<nblock; b++) {
        b_isize = blocksize[b].first;
        b_jsize = blocksize[b].second;
        fBlockPos[b+1] = fBlockPos[b]+b_isize*b_isize;
        ndata += b_isize*b_jsize;
        nr += b_isize;
        nc += b_jsize;
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
    this->fRow = nr;
    this->fCol = nc;
}

template<class TVar>
void TRMIrregularBlockDiagonal<TVar>::BuildFromMatrix(TPZMatrix<TVar> &mat) {
    if(mat.Rows() != this->Rows()) {
        cout << "TPZBlockDiagonal::BuildFromMatrix wrong data structure\n";
        return;
    }
    long nblock = fBlockSize.NElements();
    long b,nr=0,nc=0;
    for(b=0; b<nblock; b++) {
        int r,c;
        long b_isize = fBlockSize[b].first;
        long b_jsize = fBlockSize[b].second;
        long pos = fBlockPos[b];
        for(r=0; r<b_isize; r++){
            for(c=0; c<b_jsize; c++) {
                fStorage[pos+r+b_isize*c] = mat.GetVal(nr+r,nc+c);
            }
        }
        nr += b_isize;
        nc += b_jsize;
    }
}

/*******************/
/*** Constructor ***/

template<class TVar>
TRMIrregularBlockDiagonal<TVar>::TRMIrregularBlockDiagonal()
: TPZMatrix<TVar>(), fStorage(), fBlockPos(1,0),fBlockSize()
{
    
}

template<class TVar>
TRMIrregularBlockDiagonal<TVar>::TRMIrregularBlockDiagonal(TPZVec< std::pair<long, long> > &blocksize)
: TPZMatrix<TVar>(), fStorage(), fBlockPos(1,0),fBlockSize()
{
    Initialize(blocksize);
}

/********************/
/*** Constructors ***/

template<class TVar>
TRMIrregularBlockDiagonal<TVar>::TRMIrregularBlockDiagonal(TPZVec< std::pair<long, long> > &blocksizes, const TPZFMatrix<TVar> &glob)
: TPZMatrix<TVar>(), fBlockSize(blocksizes)
{
    long nblock = blocksizes.NElements();
    fBlockPos.Resize(nblock+1,0);
    long b;
    long ndata = 0;
    long nr = 0, nc = 0;
    long b_isize;
    long b_jsize;
    for(b=0; b<nblock; b++) {
        b_isize = blocksizes[b].first;
        b_jsize = blocksizes[b].second;
        fBlockPos[b+1] = fBlockPos[b]+b_isize*b_jsize;
        ndata += b_isize*b_jsize;
        nr += b_isize;
    }
    fStorage.Resize(ndata,0.);
    this->fRow = nr;
    this->fCol = nc;
    long pos, r, c;
    nr = 0;
    nc = 0;
    for(b=0; b<nblock; b++) {
        b_isize = fBlockSize[b].first;
        b_jsize = fBlockSize[b].second;
        pos = fBlockPos[b];
        for(r=0; r<b_isize; r++) {
            for(c=0; c<b_jsize; c++) {
                fStorage[pos+r+b_isize*c]= glob.GetVal(nr+r,nc+c);
            }
        }
        nr += b_isize;
        nc += b_jsize;
    }
    
}

/*********************************/
/*** Constructor( TRMIrregularBlockDiagonal& ) ***/
template<class TVar>
TRMIrregularBlockDiagonal<TVar>::TRMIrregularBlockDiagonal (const TRMIrregularBlockDiagonal<TVar> & A)
: TPZMatrix<TVar>( A.Dim(), A.Dim() ), fStorage(A.fStorage),
fBlockPos(A.fBlockPos), fBlockSize(A.fBlockSize)
{
}

/******************/
/*** Destructor ***/
template<class TVar>
TRMIrregularBlockDiagonal<TVar>::~TRMIrregularBlockDiagonal ()
{
}

/***********/
/*** Put ***/
template<class TVar>
int TRMIrregularBlockDiagonal<TVar>::Put(const long row,const long col,const TVar& value )
{
    //  cout << "TRMIrregularBlockDiagonal.Put should not be called\n";
    if ( (row >= Dim()) || (col >= Dim()) || row<0 || col<0)
    {
        cout << "TRMIrregularBlockDiagonal::Put: " << row << "," << col << "," << Dim();
        cout << "\n";
        return( 0 );
    }
    
    return( PutVal( row, col, value ) );
}

/***********/
/*** PutVal ***/
template<class TVar>
int TRMIrregularBlockDiagonal<TVar>::PutVal(const long row,const long col,const TVar& value )
{
    long b = 0;
    long nb = fBlockSize.NElements();
    if(nb==0) {
        cout << "TRMIrregularBlockDiagonal::PutVal called with parameters out of range\n";
        return -1;
    }
    long nc=0, nr=0;
    long b_isize = fBlockSize[b].first;
    long b_jsize = fBlockSize[b].second;
    while(nr+b_isize <= row && nc+b_jsize <= col && b < nb) {
        nr+=b_isize;
        nc+=b_jsize;
        b++;
        b_isize = fBlockSize[b].first;
        b_jsize = fBlockSize[b].second;
    }
    if(b==nb) {
        cout << "TRMIrregularBlockDiagonal::PutVal wrong data structure\n";
        return -1;
    }
    if(row >= nr+b_isize || col >= nc+b_jsize) {
        if(value != TVar(0.)) {
            cout << "TRMIrregularBlockDiagonal::PutVal, indices row col out of range\n";
            return -1;
        } else {
            return 0;
        }
    }
    fStorage[fBlockPos[b]+(row-nr)+b_isize*(col-nc)] = value;
    return 0;
}



/***********/
/*** Get ***/

template<class TVar>
const TVar&
TRMIrregularBlockDiagonal<TVar>::Get(const long row,const long col ) const
{
    if ( (row >= Dim()) || (col >= Dim()) )
        TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "TRMIrregularBlockDiagonal::Get <indices out of band matrix range>" );
    
    return( GetVal( row, col ) );
}

template<class TVar>
TVar &
TRMIrregularBlockDiagonal<TVar>::operator()(const long row, const long col) {
    
    long b = 0;
    long nb = fBlockSize.NElements();
    if(nb==0) {
        cout << "TRMIrregularBlockDiagonal::operator() called with parameters out of range\n";
        static TVar zero = 0.;
        return zero;
    }
    long nr=0, nc=0;
    long b_isize = fBlockSize[b].first;
    long b_jsize = fBlockSize[b].second;
    while(nr+b_isize <= row && nc+b_jsize <= col && b < nb) {
        nr+=b_isize;
        nc+=b_jsize;
        b++;
        b_isize = fBlockSize[b].first;
        b_jsize = fBlockSize[b].second;
    }
    if(b==nb) {
        cout << "TRMIrregularBlockDiagonal::operator() wrong data structure\n";
        static TVar zero = 0.;
        return zero;
    }
    if(row >= nr+b_isize || col >= nc+b_jsize) {
        cout << "TRMIrregularBlockDiagonal::operator(), indices row col out of range\n";
        static TVar zero = 0.;
        return zero;
    }
    return fStorage[fBlockPos[b]+(row-nr)+b_isize*(col-nc)];
}

/***********/
/*** GetVal ***/
template<class TVar>
const TVar &TRMIrregularBlockDiagonal<TVar>::GetVal(const long row,const long col ) const
{
    long b = 0;
    long nb = fBlockSize.NElements();
    if(nb==0) {
        cout << "TPZBlockDiagonal::GetVal called with parameters out of range\n";
    }
    long nr=0, nc=0;
    long b_isize = fBlockSize[b].first;
    long b_jsize = fBlockSize[b].second;
    while(nr+b_isize <= row && nc+b_jsize <= col && b < nb) {
        nr+=b_isize;
        nc+=b_jsize;
        b++;
        b_isize = fBlockSize[b].first;
        b_jsize = fBlockSize[b].second;
    }
    if(b==nb) {
        cout << "TPZBlockDiagonal::GetVal wrong data structure\n";
    }
    if(row >= nr+b_isize || col >= nc+b_jsize) {
        //cout << "TPZBlockDiagonal::GetVal, indices row col out of range\n";
        static TVar zero = 0.;
        return zero;
    }
    return fStorage[fBlockPos[b]+(row-nr)+b_isize*(col-nc)];
}


/******** Operacoes com MATRIZES GENERICAS ********/

/*******************/
/*** MultiplyAdd ***/
//
//  perform a multiply add operation to be used by iterative solvers
//

template<class TVar>
void TRMIrregularBlockDiagonal<TVar>::MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
                                     const TVar alpha,const TVar beta ,const int opt) const {
    // Computes z = beta * y + alpha * opt(this)*x
    //          z and x cannot overlap in memory
    
    if ((!opt && this->Cols() != x.Rows()) || this->Rows() != x.Rows())
        TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "TRMIrregularBlockDiagonal::MultAdd <matrixs with incompatible dimensions>" );
    if(x.Cols() != y.Cols() || x.Cols() != z.Cols() || x.Rows() != y.Rows() || x.Rows() != z.Rows()) {
        TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"TRMIrregularBlockDiagonal::MultAdd incompatible dimensions\n");
    }
    
    this->PrepareZ(y,z,beta,opt);
    long xcols = x.Cols();
    long nb= fBlockSize.NElements();
    long ic, b, nr=0, nc=0;
    long b_isize, b_jsize, r, c;
    if(opt == 0) {
        for (ic = 0; ic < xcols; ic++) {
			nr=0, nc=0;
            for(b=0; b<nb; b++) {
                b_isize = fBlockSize[b].first;
                b_jsize = fBlockSize[b].second;
                long pos = fBlockPos[b];
                for(r=0; r<b_isize; r++) {
                    for(c=0; c<b_jsize; c++) {
                        z(nr+r,ic) += alpha*fStorage[pos+r+b_isize*c]*x.GetVal((nc+c),ic);
                    }
                }
                nr += b_isize;
                nc += b_jsize;
            }
        }
    } else {
        cout << "xcols \t" << xcols << "\n";
        for (ic = 0; ic < xcols; ic++) {
			nr=0, nc=0;
            for(b=0; b<nb; b++) {
                b_isize = fBlockSize[b].first;
                b_jsize = fBlockSize[b].second;
                long pos = fBlockPos[b];
                for(r=0; r<b_isize; r++) {
                    for(c=0; c<b_jsize; c++) {
                        z(nr+r,ic) += alpha*fStorage[pos+r+b_isize*c]*x.GetVal((nc+c),ic);
                    }
                }
                nr += b_isize;
                nc += b_jsize;
            }
        }
    }
}


/***************/
/**** Zero ****/
template<class TVar>
int TRMIrregularBlockDiagonal<TVar>::Zero()
{
    
    fStorage.Fill(0.,0);
    this->fDecomposed = 0;
    
    return( 1 );
}



/********************/
/*** Transpose () ***/
template<class TVar>
void TRMIrregularBlockDiagonal<TVar>::Transpose (TPZMatrix<TVar> *const T) const
{
    T->Resize( Dim(), Dim() );
    
    long b, nr=0, nc=0, pos;
    long b_isize, b_jsize, r, c;
    long nb = fBlockSize.NElements();
    for ( b=0; b<nb; b++) {
        pos= fBlockPos[b];
        b_isize = fBlockSize[b].first;
        b_jsize = fBlockSize[b].second;
        for(r=0; r<b_isize; r++) {
            for(c=0; c<b_jsize; c++) {
                T->PutVal(nr+r,nc+c,fStorage[pos+c+r*b_jsize]);
            }
        }
        nr += b_isize;
        nc += b_jsize;
    }
}


/*****************/
/*** Decompose_LU ***/
//fElem[ fBand * (2*row + 1) + col ]
template<class TVar>
int TRMIrregularBlockDiagonal<TVar>::Decompose_LU(std::list<long> &singular)
{
    return Decompose_LU(singular);
}

template<class TVar>
int TRMIrregularBlockDiagonal<TVar>::Decompose_LU()
{
    
    LOGPZ_DEBUG(logger, "TRMIrregularBlockDiagonal::Decompose_LU");
    
    if (  this->fDecomposed && this->fDecomposed == ELU) {
        return ELU;
    } else if(this->fDecomposed) {
        TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"TRMIrregularBlockDiagonal::Decompose_LU is already decomposed with other scheme");
    }
    
    long b,nb,pos;
	long b_isize, b_jsize;
    nb = fBlockSize.NElements();
    for(b=0;b<nb; b++) {
        
        pos = fBlockPos[b];
        b_isize = fBlockSize[b].first;
        b_jsize = fBlockSize[b].second;
        if(!b_isize || !b_jsize) continue;
        if (b_isize != b_jsize) {
            TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"TRMIrregularBlockDiagonal::Decompose_LU, rectangular block matrix, set all blocks square");
            DebugStop();
        }
        
#ifdef LOG4CXX
        if (logger->isDebugEnabled())
        {
            std::stringstream mess;
            mess << "TRMIrregularBlockDiagonal::Decompose_LU() - bsize = " << bsize << ", bsize*bsize = " << bsize*bsize;
            LOGPZ_DEBUG(logger,mess.str());
        }
#endif
        
        TPZFMatrix<TVar> temp(b_isize,b_jsize,&fStorage[pos],b_isize*b_jsize);
        std::list<long> singular;
        temp.Decompose_LU(singular);
    }
    this->fDecomposed = ELU;
    return 1;
}

template<class TVar>
int
TRMIrregularBlockDiagonal<TVar>::Substitution( TPZFMatrix<TVar> *B) const
{
    if(this->fDecomposed != ELU) {
        TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"TRMIrregularBlockDiagonal::Decompose_LU is decomposed with other scheme");
    }
    
    long b,nb,pos,b_isize,b_jsize,eq;
    nb = fBlockSize.NElements();
    long c, nc = B->Cols();
    for(c=0; c<nc; c++) {
        eq = 0;
        for(b=0;b<nb; b++) {
            pos = fBlockPos[b];
            b_isize = fBlockSize[b].first;
            b_jsize = fBlockSize[b].second;
            if(!b_isize || !b_jsize) continue;
            if (b_isize != b_jsize) {
                TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"TPZBlockDiagonal::Substitution, rectangular block matrix, set all blocks square");
                DebugStop();
            }
            TPZFMatrix<TVar> BTemp(b_jsize,1,&(B->operator()(eq,c)),b_jsize);
            TVar *ptr = fStorage.begin()+pos;
            TPZFMatrix<TVar>::Substitution(ptr,b_jsize,&BTemp);
            eq+= b_jsize;
        }
    }
    return 1;
}

/************************** Private **************************/

/*************/
/*** Clear ***/
template<class TVar>
int TRMIrregularBlockDiagonal<TVar>::Clear()
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
int TRMIrregularBlockDiagonal<TVar>::main() {
    
    cout << "Entering the main program\n";
    cout.flush();
    TPZFMatrix<TVar> ref(7,7,0.);
    int r,c;
    for(r=0; r<7; r++) {
        for(c=0; c<7; c++) {
            ref(r,c) = ((TVar)(float)(5+r*c+3*r));
        }
        ref(r,r) += (TVar)1000;
    }
    TPZVec< pair<long, long> > blocksize(3);
    blocksize[0].first = 2;
    blocksize[0].second = 2;
    blocksize[1].first = 4;
    blocksize[1].second = 4;
    blocksize[2].first = 1;
    blocksize[2].second = 1;
    TRMIrregularBlockDiagonal bd1(blocksize,ref);
    TRMIrregularBlockDiagonal bd2(bd1);
    ref.Print("original matrix",std::cout);
    bd1.Solve_LU(&ref);
    bd1.Solve_LU(&ref);
    ref.Print("after inverting the diagonal",std::cout);
    TPZFMatrix<TVar> ref2;
    bd2.Multiply(ref,ref2);
    bd2.Multiply(ref2,ref);
    ref.Print("restoring the original matrix",std::cout);
    return 1;
    
}

template<class TVar>
void TRMIrregularBlockDiagonal<TVar>::Print(const char *msg, std::ostream &out, const MatrixOutputFormat format) const {
    
    if(format != EFormatted)
    {
        TPZMatrix<TVar>::Print(msg,out,format);
        return;
    }
    out << "TRMIrregularBlockDiagonal matrix ";
    if(msg) out << msg;
    out  << std::endl;
    
    long nblock = fBlockSize.NElements();
    out << "Number of blocks " << nblock << std::endl;
	long b,b_isize,b_jsize,pos;
    for(b=0; b<nblock; b++) {
        b_isize = fBlockSize[b].first;
        b_jsize = fBlockSize[b].second;
        out << "block number " << b << " size : " << b_isize << " x " << b_jsize << std::endl;
        long r,c;
        pos = fBlockPos[b];
        for(c=0; c<b_isize; c++) {
            for(r=0; r<b_jsize ; r++) {
                out << fStorage[pos+r+b_isize*c] << ' ';
            }
            out << std::endl;
        }
    }
}

/**
 * Updates the values of the matrix based on the values of the matrix
 */
template<class TVar>
void TRMIrregularBlockDiagonal<TVar>::UpdateFrom(TPZAutoPointer<TPZMatrix<TVar> > mat)
{
    if(!mat)
    {
        cout << "TRMIrregularBlockDiagonal::UpdateFrom" << " called with zero argument\n";
        return;
    }
    this->fDecomposed = ENoDecompose;
    long nblock = fBlockSize.NElements();
    long b,b_isize,b_jsize,pos,firstrow = 0, firstcol = 0;
    for(b=0; b<nblock; b++) {
        b_isize = fBlockSize[b].first;
        b_jsize = fBlockSize[b].second;
        //    int r,c;
        pos = fBlockPos[b];
        TPZFMatrix<TVar> block(b_isize,b_jsize,&fStorage[pos],b_isize*b_jsize);
        mat->GetSub(firstrow,firstcol,b_isize,b_jsize,block);
        firstrow += b_isize;
        firstcol += b_isize;
    }
}

/** Fill the matrix with random values (non singular matrix) */
template<class TVar>
void TRMIrregularBlockDiagonal<TVar>::AutoFill(long neq, long jeq, int symmetric) {
    
    if (neq != jeq) {
        DebugStop();
    }
    TPZStack<std::pair<long,long> > blsizes;
    long totalsize = 0;
    while (totalsize < neq) {
        long blsize = (neq*rand())/RAND_MAX;
        blsize = blsize < neq-totalsize ? blsize : neq-totalsize;
        blsizes.Push(std::make_pair(blsize,blsize));
        totalsize += blsize;
    }
    Initialize(blsizes);
    // Initialize the blocksizes!!
    long b, bsize, eq = 0, pos;
    long nb = fBlockSize.NElements(), r, c;
    for ( b=0; b<nb; b++) {
        pos= fBlockPos[b];
        bsize = fBlockSize[b].first;
        for(c=0; c<bsize; c++) {
            float sum = 0.;
            r=0;
            if (symmetric == 1) {
                for (r=0; r<c; r++) {
                    fStorage[pos+c+r*bsize]=fStorage[pos+r+c*bsize];
                    sum += fabs(fStorage[pos+r+c*bsize]);
                }
            }
            for(; r<bsize; r++) {
                float val = ((float)rand())/RAND_MAX;
                fStorage[pos+c+r*bsize] = (TVar)(val);
                if(c!= r) sum += fabs(val);
            }
            if (fabs(fStorage[pos+c+c*bsize]) < sum) {
                fStorage[pos+c+c*bsize] = (TVar)(sum + (float)1.);
            }
        }
        eq += bsize;
    }
}

template class TRMIrregularBlockDiagonal<float>;
template class TRMIrregularBlockDiagonal<double>;
template class TRMIrregularBlockDiagonal<long double>;

template class TRMIrregularBlockDiagonal<std::complex<float> >;
template class TRMIrregularBlockDiagonal<std::complex<double> >;
template class TRMIrregularBlockDiagonal<std::complex<long double> >;