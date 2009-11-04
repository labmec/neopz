//
//
// Author: MISAEL LUIS SANTANA MANDUJANO.
//
// File:   tmatred.c
//
// Class:  TPZMatRed
//
// Obs.: SubestruturaÔøΩo simples de um sistema de equacoes.
//       So'para matrizes quadradas
// Versao: 04 / 1996.
//

#include <stdlib.h>
#include <stdio.h>
#include <fstream>
using namespace std;


#include "pzmatred.h"
#include "pzfmatrix.h"

#include <sstream>
#include "pzlog.h"
#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.matrix.tpzmatred"));
#endif

//REAL TPZMatrix::gZero = 0.;

//ofstream out("saida");

/*************************** Public ***************************/

/******************/
/*** Construtor ***/
template<class TSideMatrix>
TPZMatRed<TSideMatrix>::TPZMatRed () : TPZMatrix( 0, 0 ), fK01(0,0),fK10(0,0),fK11(0,0),fF0(0,0),fF1(0,0)
{
  fDim0=0;
  fDim1=0;
  fDecomposeType=ENoDecompose;
  fF0IsComputed=0;
  fK11IsReduced=0;
	fK01IsComputed = 0;
  fF1IsReduced=0;

}


template<class TSideMatrix>
TPZMatRed<TSideMatrix>::TPZMatRed( int dim, int dim00 ):TPZMatrix( dim,dim ), fK01(dim00,dim-dim00,0.), 
	fK10(dim-dim00,dim00,0.), fK11(dim-dim00,dim-dim00,0.), fF0(dim00,1,0.),fF1(dim-dim00,1,0.) 
{

  if(dim<dim00) TPZMatrix::Error(__PRETTY_FUNCTION__,"dim k00> dim");
  fDim0=dim00;
  fDim1=dim-dim00;
  fDecomposeType=ENoDecompose;
  fF0IsComputed=0;
  fK11IsReduced=0;
	fK01IsComputed = 0;
  fF1IsReduced=0;

}


template<class TSideMatrix>
TPZMatRed<TSideMatrix>::TPZMatRed(const TPZMatRed &cp) : TPZMatrix(cp), fK01(cp.fK01), fK10(cp.fK10), fK11(cp.fK11), fF0(cp.fF0), fF1(cp.fF1)
{
  fDim0=cp.fDim0;
  fDim1=cp.fDim1;
  fDecomposeType=cp.fDecomposeType;
  fF0IsComputed=cp.fF0IsComputed;
  fK11IsReduced=cp.fK11IsReduced;
	fK01IsComputed = cp.fK01IsComputed;
  fF1IsReduced=fF1IsReduced;

  if(cp.fK00) fK00 = cp.fK00->Clone();
}


template<class TSideMatrix>
TPZMatRed<TSideMatrix>::~TPZMatRed(){
}


template<class TSideMatrix>
int TPZMatRed<TSideMatrix>::IsSimetric() const {
  if(fK00) return fK00->IsSimetric();
  return 0;
}


template<class TSideMatrix>
void TPZMatRed<TSideMatrix>::Simetrize() {
  // considering fK00 is simetric, only half of the object is assembled.
  // this method simetrizes the matrix object

  if(!fK00 || !fK00->IsSimetric()) return;
  fK01.Transpose(&fK10);
  int row,col;
  for(row=0; row<fDim1; row++) {
    for(col=row+1; col<fDim1; col++) {
      (fK11)(col,row) = (fK11)(row,col);
    }
  }
}

template<class TSideMatrix>
int
TPZMatRed<TSideMatrix>::PutVal(const int r,const int c,const REAL& value ){
  int row(r),col(c);
  if (IsSimetric() && row > col ) Swap( &row, &col );
  if (row<fDim0 &&  col<fDim0)  fK00->PutVal(row,col,value);
  if (row<fDim0 &&  col>=fDim0)  fK01.PutVal(row,col-fDim0,value);
  if (row>=fDim0 &&  col<fDim0)  fK10.PutVal(row-fDim0,col,value);
  if (row>=fDim0 &&  col>=fDim0)  fK11.PutVal(row-fDim0,col-fDim0,value);


  return( 1 );
}


template<class TSideMatrix>
const REAL&
TPZMatRed<TSideMatrix>::GetVal(const int r,const int c ) const {
  int row(r),col(c);

  if (IsSimetric() && row > col ) Swap( &row, &col );
  if (row<fDim0 &&  col<fDim0)  return ( fK00->GetVal(row,col) );
  if (row<fDim0 &&  col>=fDim0)  return ( fK01.GetVal(row,col-fDim0) );
  if (row>=fDim0 &&  col<fDim0)  return ( fK10.GetVal(row-fDim0,col) );
  return (fK11.GetVal(row-fDim0,col-fDim0) );

}

template<class TSideMatrix>
REAL&
TPZMatRed<TSideMatrix>::s(const int r,const int c ) {
	int row(r),col(c);
	
	if (r < fDim0 && IsSimetric() && row > col ) Swap( &row, &col );
	if (row<fDim0 &&  col<fDim0)  return ( fK00->s(row,col) );
	if (row<fDim0 &&  col>=fDim0)  return ( fK01.s(row,col-fDim0) );
	if (row>=fDim0 &&  col<fDim0)  return ( fK10.s(row-fDim0,col) );
	return (fK11.s(row-fDim0,col-fDim0) );
	
}


template<class TSideMatrix>
void
TPZMatRed<TSideMatrix>::SetSolver(TPZAutoPointer<TPZMatrixSolver> solver)
{
	fK00=solver->Matrix();
	fSolver = solver;
}

template<class TSideMatrix>
void
TPZMatRed<TSideMatrix>::SetK00(TPZAutoPointer<TPZMatrix> K00)
{
	fK00=K00;
}

template<class TSideMatrix>
void
TPZMatRed<TSideMatrix>::SetF(const TPZFMatrix & F)
{

  int FCols=F.Cols(),c,r,r1;

  fF0.Redim(fDim0,FCols);
  fF1.Redim(fDim1,FCols);

  for(c=0; c<FCols; c++){
    r1=0;
    for(r=0; r<fDim0; r++){
      fF0.PutVal( r,c,F.GetVal(r,c) ) ;
    }
    //aqui r=fDim0
    for( ;r<Rows(); r++){
      fF1.PutVal( r1++,c,F.GetVal(r,c) );
    }
  }
  fF1IsReduced = 0;
  fF0IsComputed = 0;
}

template<class TSideMatrix>
const TPZFMatrix&
TPZMatRed<TSideMatrix>::F1Red()
{
	if (!fDim0 || fF1IsReduced)  return (fF1);
	if(! fF0IsComputed)
	{
		fSolver->Solve(fF0,fF0);
		fF0IsComputed = 1;
	}
	
	//make [F1]=[F1]-[K10][K0]
	fK10.MultAdd((fF0),(fF1),(fF1),-1,1);
	fF1IsReduced=1;
	return (fF1);
}

template<class TSideMatrix>
const TPZFMatrix&
TPZMatRed<TSideMatrix>::K11Red()
{
	if (!fDim0 || fK11IsReduced)  {
		//Simetrize();
		return (fK11);
	}
	
	if(!fK01IsComputed)
	{
		Simetrize();
		fSolver->Solve(fK01,fK01);
		fK01IsComputed = 1;
	}

	//make [K11]=[k11]-[k10][k01]
	fK10.MultAdd(fK01,(fK11),(fK11),-1,1);
	fK11IsReduced=1;
	return (fK11);


}

#include "tpzverysparsematrix.h"

template<>
const TPZFMatrix&
TPZMatRed<TPZVerySparseMatrix>::K11Red()
{
	std::cout << __PRETTY_FUNCTION__ << " should never be called\n";	
	static TPZFMatrix temp;
	return temp;
}


template<class TSideMatrix>
void
TPZMatRed<TSideMatrix>::U1(TPZFMatrix & F)
{

  K11Red();
  F1Red();
  F=(fF1);
  fK11.SolveDirect( F ,ELU);


}

template<>
void
TPZMatRed<TPZVerySparseMatrix>::UGlobal(const TPZFMatrix & U1, TPZFMatrix & result)
{
	//[u0]=[A00^-1][F0]-[A00^-1][A01]
    if( !fF0IsComputed ){
		//compute [F0]=[A00^-1][F0]
		fSolver->Solve(fF0,fF0);
		fF0IsComputed=1;
	}
	
	if(!fK01IsComputed)
	{
		TPZFMatrix k01(fK01);
		fSolver->Solve(k01,k01);
		fK01 = k01;
		fK01IsComputed = 1;
	}
	
	//make [u0]=[F0]-[U1]
	TPZFMatrix u0( fF0.Rows() , fF0.Cols() );
	fK01.MultAdd(U1,(fF0),u0,-1,0);
	//	fSolver->Solve(u0,u0);
	//	u0 += fF0;
	//compute result
	
	result.Redim( Rows(),fF0.Cols() );
	int c,r,r1;
	
	for(c=0; c<fF0.Cols(); c++)
	{
		r1=0;
		for(r=0; r<fDim0; r++)
		{
			result.PutVal( r,c,u0.GetVal(r,c) ) ;
		}
		//aqui r=fDim0
		for( ;r<Rows(); r++)
		{
			result.PutVal( r,c,U1.GetVal(r1++,c) );
		}
	}
}

template<class TSideMatrix>
void
TPZMatRed<TSideMatrix>::UGlobal(const TPZFMatrix & U1, TPZFMatrix & result)
{
  //[u0]=[A00^-1][F0]-[A00^-1][A01]
    if( !fF0IsComputed ){
		//compute [F0]=[A00^-1][F0]
		fSolver->Solve(fF0,fF0);
		fF0IsComputed=1;
	}
	
	if(!fK01IsComputed)
	{
		fSolver->Solve(fK01,fK01);
		fK01IsComputed = 1;
	}

	//make [u0]=[F0]-[U1]
	TPZFMatrix u0( fF0.Rows() , fF0.Cols() );
	fK01.MultAdd(U1,(fF0),u0,-1,1);
//	fSolver->Solve(u0,u0);
//	u0 += fF0;
	//compute result
#ifdef LOG4CXX
	if(logger->isDebugEnabled())
	{
		std::stringstream sout;
		fF0.Print("fF0 ",sout);
		u0.Print("u0 " ,sout);
		LOGPZ_DEBUG(logger,sout.str())   

	}
#endif

	result.Redim( Rows(),fF0.Cols() );
	int c,r,r1;

	for(c=0; c<fF0.Cols(); c++)
	{
		r1=0;
		for(r=0; r<fDim0; r++)
		{
			result.PutVal( r,c,u0.GetVal(r,c) ) ;
		}
		//aqui r=fDim0
		for( ;r<Rows(); r++)
		{
			result.PutVal( r,c,U1.GetVal(r1++,c) );
		}
	}
}


template<class TSideMatrix>
void
TPZMatRed<TSideMatrix>::Print(const char *name , std::ostream &out ,const MatrixOutputFormat form ) const
{


  if(form != EInputFormat) {
    out << "Writing matrix 'TPZMatRed<TSideMatrix>::" << name;
    out << "' (" << Dim() << " x " << Dim() << "):\n";

    fK00->Print("K(0,0)",out,form);
    fK01.Print("K(0,1)",out,form);
    fK10.Print("K(1,0)",out,form);
    fK11.Print("K(1,1)",out,form);


      fF0.Print("F(0)",out,form);
      fF1.Print("F(1)",out,form);

    out << "\n\n";
  } else {
    TPZMatrix::Print(name,out,form);
  }

}

template<class TSideMatrix>
int TPZMatRed<TSideMatrix>::Substitution(TPZFMatrix *B) const{
	switch(fDecomposeType) {
		case ENoDecompose:
			default:
				TPZMatrix::Error(__PRETTY_FUNCTION__, "TPZMatRed<TSideMatrix>::Substitution called without initialized solver\n");
				return 0;

		case ELU:
			return( fK00->Substitution(B) );

		case ECholesky:
			return ( fK00->Subst_Forward( B ) && fK00->Subst_Backward( B )  );

		case ELDLt:
			return ( fK00->Subst_LForward( B )
				 && fK00->Subst_Diag( B ) && fK00->Subst_LBackward( B ) );
	}
}

template<class TSideMatrix>
int TPZMatRed<TSideMatrix>::Redim(int dim, int dim00){
	if(dim<dim00) TPZMatrix::Error(__PRETTY_FUNCTION__,"dim k00> dim");
	if(fK00) fK00->Redim(dim00,dim00);
	if(fF0)  delete fF0;
	if(fF1)  delete fF1;

	fDim0=dim00;
	fDim1=dim-dim00;
	fDecomposeType=ENoDecompose;
	fF0IsComputed=0;
	fK11IsReduced=0;
	fF1IsReduced=0;

	fK01.Redim(fDim0,fDim1);
	fK10.Redim(fDim1,fDim0);
	fK11.Redim(fDim1,fDim1);

	fF0=NULL;
	fF1=NULL;
	fRow = dim;
	fCol = dim;
	return 0;
}


template<class TSideMatrix>
int TPZMatRed<TSideMatrix>::Zero(){
	if(fK00) fK00->Zero();
	fK01.Zero();
	fK10.Zero();
	fK11.Zero();
	fF0.Zero();
	fF1.Zero();
	fF0IsComputed=0;
	fK11IsReduced=0;
	fF1IsReduced=0;
	return 0;
}


template<class TSideMatrix>
void TPZMatRed<TSideMatrix>::MultAdd(const TPZFMatrix &x,
                        const TPZFMatrix &y, TPZFMatrix &z,
                        const REAL alpha,const REAL beta,
                        const int opt,const int stride) const
{
#warning Not functional yet. Still need to Identify all the variables
//   int i = 0;
//   int j = 0;

	if(!fIsReduced)
	{
		TPZMatrix::MultAdd(x,y,z,alpha,beta,opt,stride);
		return;
	}
  /**
   * It computes z = beta * y + alpha * opt(this)*x but z and x can not overlap in memory.
   * @param x Is x on the above operation
   * @param y Is y on the above operation
   * @param z Is z on the above operation
   * @param alpha Is alpha on the above operation
   * @param beta Is beta on the above operation
   * @param opt Indicates if is Transpose or not
   * @param stride Indicates n/N where n is dimension of the right hand side vector and N is matrix dimension
   */
	this->PrepareZ(y,z,beta,opt,stride);

	if(!opt)
	{
		TPZFMatrix l_Res(fK01.Rows(), x.Cols(), 0);
		fK01.Multiply(x,l_Res,0,1);
		fSolver->Solve(l_Res,l_Res);
#ifdef LOG4CXX
		if(logger->isDebugEnabled())
		{
			std::stringstream sout;
			l_Res.Print("Internal solution",sout);
			LOGPZ_DEBUG(logger,sout.str())
		}
#endif
		TPZFMatrix l_ResFinal(fK11.Rows(), x.Cols(), 0);
		fK11.Multiply(x,l_ResFinal,0,1);
#ifdef LOG4CXX
		if(logger->isDebugEnabled())
		{
			std::stringstream sout;
			l_ResFinal.Print("Intermediate product l_ResFinal",sout);
			LOGPZ_DEBUG(logger,sout.str())
		}
#endif
		fK10.MultAdd(l_Res,l_ResFinal,z,-alpha,alpha,opt,stride);
#ifdef LOG4CXX
		if(logger->isDebugEnabled())
		{
			std::stringstream sout;
			z.Print("Final result z ",sout);
			LOGPZ_DEBUG(logger,sout.str())
		}
#endif
	}
	else
	{
		DebugStop();
	}
}

template<class TSideMatrix>
void TPZMatRed<TSideMatrix>::Write(TPZStream &buf, int withclassid)
{
  TPZMatrix::Write(buf, withclassid);
  {//Ints
    buf.Write(&this->fDim0, 1);
    buf.Write(&this->fDim1, 1);
    int dcType = (int)this->fDecomposeType;
    buf.Write(&dcType, 1);
  }
  {//chars
    buf.Write(&this->fDecomposed, 1);
    buf.Write(&this->fDefPositive, 1);
    buf.Write(&this->fF0IsComputed, 1);
    buf.Write(&this->fF1IsReduced, 1);
    buf.Write(&this->fIsReduced, 1);
    buf.Write(&this->fK01IsComputed, 1);
    buf.Write(&this->fK11IsReduced, 1);
  }
  {//Aggregates
    this->fF0.Write(buf, 0);
    this->fF1.Write(buf, 0);
    if(this->fK00)
    {
      this->fK00->Write(buf, 0);
    }
    else
    {
      int flag = -1;
      buf.Write(&flag, 1);
    }
    this->fK01.Write(buf, 0);
    this->fK10.Write(buf, 0);
    this->fK11.Write(buf, 0);
    if(fSolver)
    {
      if(fSolver->Matrix() != fK00)
      {
        std::cout << "Error\n";
      }
      else
      {
        fSolver->Write(buf, 0);
        //TODO Enviar o solver, atenção com a Matrix do Solver;
      }

    }
    else
    {
      int flag = -1;
      buf.Write(&flag, 1);
    }

  }

}
template<class TSideMatrix>
void TPZMatRed<TSideMatrix>::Read(TPZStream &buf, void *context)
{
  TPZMatrix::Read(buf, context);
  {//Ints
    buf.Read(&this->fDim0, 1);
    buf.Read(&this->fDim1, 1);
    int dcType = 0;
    buf.Read(&dcType, 1);
    this->fDecomposeType = (DecomposeType) dcType;
  }
  {//chars
    buf.Read(&this->fDecomposed, 1);
    buf.Read(&this->fDefPositive, 1);
    buf.Read(&this->fF0IsComputed, 1);
    buf.Read(&this->fF1IsReduced, 1);
    buf.Read(&this->fIsReduced, 1);
    buf.Read(&this->fK01IsComputed, 1);
    buf.Read(&this->fK11IsReduced, 1);
  }
  {//Aggregates
    this->fF0.Read(buf, 0);
    this->fF1.Read(buf, 0);
    if(!fK00)
    {
      TSideMatrix * side = new TSideMatrix;
      side->Read(buf, 0);
      fK00 = side;
    }
    this->fK01.Read(buf, 0);
    this->fK10.Read(buf, 0);
    this->fK11.Read(buf, 0);
//    if(!fSolver)
//    {
//      TPZMatrixSolver * solver = new TPZMatrixSolver;
//      solver->Read(buf, 0);
//      fSolver = solver;
//    }
  }
}

/*************************** Private ***************************/

/*************/
/*** Error ***/

/*int
TPZMatRed<TSideMatrix>::Error(const char *msg ,const char *msg2)
{
  ostringstream out;
  out << "TPZMatRed<TSideMatrix>::" << msg << msg2 << ".\n";
  LOGPZ_ERROR (logger, out.str().c_str());
  exit( 1 );
  return 0;
}*/



#include "tpzverysparsematrix.h"

template class TPZMatRed<TPZVerySparseMatrix>;
template class TPZMatRed<TPZFMatrix>;

template <>
int TPZMatRed<TPZVerySparseMatrix>::ClassId()
{
  return TPZMATRED_VERYSPARSE_ID;
}
template <>
int TPZMatRed<TPZFMatrix>::ClassId()
{
  return TPZMATRED_FMATRIX_ID;
}

template class TPZRestoreClass<TPZMatRed<TPZVerySparseMatrix>, TPZMATRED_VERYSPARSE_ID>;
template class TPZRestoreClass<TPZMatRed<TPZFMatrix>, TPZMATRED_FMATRIX_ID>;
