//
//
// Author: MISAEL LUIS SANTANA MANDUJANO.
//
// File:   tmatred.c
//
// Class:  TPZMatRed
//
// Obs.: Subestrutura€ao simples de um sistema de equacoes.
//       So'para matrizes quadradas
// Versao: 04 / 1996.
//

#include <stdlib.h>
#include <stdio.h>
#include <fstream>
using namespace std;


#include "pzmatred.h"
#include "pzfmatrix.h"

//REAL TPZMatrix::gZero = 0.;

//ofstream out("saida");

/*************************** Public ***************************/

/******************/
/*** Construtor ***/

TPZMatRed::TPZMatRed () : TPZMatrix( 0, 0 )
{
  fK00=NULL;
  fK01=NULL;
  fK10=NULL;
  fK11=NULL;
  fF0 =NULL;
  fF1 =NULL;
  fDim0=0;
  fDim1=0;
  fDecomposeType=ENoDecompose;
  fF0IsComputed=0;
  fK01IsComputed=0;
  fK11IsReduced=0;
  fF1IsReduced=0;

}



TPZMatRed::TPZMatRed( int dim, int dim00 ):TPZMatrix( dim,dim ){

  if(dim<dim00) Error("dim k00> dim");
  fDim0=dim00;
  fDim1=dim-dim00;
  fDecomposeType=ENoDecompose;
  fF0IsComputed=0;
  fK01IsComputed=0;
  fK11IsReduced=0;
  fF1IsReduced=0;

  fK00= 0; //new TPZFMatrix(fDim0,fDim0,0.);
  fK01=new TPZFMatrix(fDim0,fDim1,0.);
  fK10=new TPZFMatrix(fDim1,fDim0,0.);
  fK11=new TPZFMatrix(fDim1,fDim1,0.);


  fF0=NULL;
  fF1=NULL;

}

TPZMatRed::~TPZMatRed(){
	if(fK00) delete fK00;
	if(fK01) delete fK01;
	if(fK10) delete fK10;
	if(fK11) delete fK11;
	if(fF0)  delete fF0;
	if(fF1)  delete fF1;
}

int TPZMatRed::IsSimetric() const {
  if(fK00) return fK00->IsSimetric();
  return 0;
}


void TPZMatRed::Simetrize() {
  // considering fK00 is simetric, only half of the object is assembled.
  // this method simetrizes the matrix object

  if(!fK00 || !fK00->IsSimetric()) return;
  fK01->Transpose(fK10);
  int row,col;
  for(row=0; row<fDim1; row++) {
    for(col=row+1; col<fDim1; col++) {
      (*fK11)(col,row) = (*fK11)(row,col);
    }
  }
}

int
TPZMatRed::PutVal(const int r,const int c,const REAL& value ){
  int row(r),col(c);
  if (IsSimetric() && row > col ) Swap( &row, &col );
  if (row<fDim0 &&  col<fDim0)  fK00->PutVal(row,col,value);
  if (row<fDim0 &&  col>=fDim0)  fK01->PutVal(row,col-fDim0,value);
  if (row>=fDim0 &&  col<fDim0)  fK10->PutVal(row-fDim0,col,value);
  if (row>=fDim0 &&  col>=fDim0)  fK11->PutVal(row-fDim0,col-fDim0,value);


  return( 1 );
}
const REAL&
TPZMatRed::GetVal(const int r,const int c ) const {
  int row(r),col(c);

  if (IsSimetric() && row > col ) Swap( &row, &col );
  if (row<fDim0 &&  col<fDim0)  return ( fK00->GetVal(row,col) );
  if (row<fDim0 &&  col>=fDim0)  return ( fK01->GetVal(row,col-fDim0) );
  if (row>=fDim0 &&  col<fDim0)  return ( fK10->GetVal(row-fDim0,col) );
  return (fK11->GetVal(row-fDim0,col-fDim0) );

}

void
TPZMatRed::SetK00(TPZMatrix *const k00)
{
  if(fK00) delete fK00;
  k00->Resize(fDim0,fDim0);
  fK00=k00;
}

void
TPZMatRed::SetF(const TPZFMatrix & F)
{

  int FCols=F.Cols(),c,r,r1;

  if(fF0) delete fF0;
  fF0=NULL;
  if(fF1) delete fF1;
  fF1=NULL;

  fF0=new TPZFMatrix(fDim0,FCols,0.);
  fF1=new TPZFMatrix(fDim1,FCols,0.);

  for(c=0; c<FCols; c++){
    r1=0;
    for(r=0; r<fDim0; r++){
      fF0->PutVal( r,c,F.GetVal(r,c) ) ;
    }
    //aqui r=fDim0
    for( ;r<Rows(); r++){
      fF1->PutVal( r1++,c,F.GetVal(r,c) );
    }


  }
  fF1IsReduced = 0;
  fF0IsComputed = 0;
}
const TPZFMatrix&
TPZMatRed::F1Red()
{
  if (!fDim0 || fF1IsReduced)  return (*fF1);

  if (fK00->IsDecomposed()) {
    if( !fF0IsComputed ){
      Substitution((fF0));
      fF0IsComputed=1;
    }
  }
  else{
    Simetrize();
    fK00->SolveDirect((*fF0),fDecomposeType);
    fF0IsComputed=1;
  }


  //make [F1]=[F1]-[K10][K0]
  fK10->MultAdd((*fF0),(*fF1),(*fF1),-1,1);
  fF1IsReduced=1;
  return (*fF1);



}

const TPZFMatrix&
TPZMatRed::K11Red()
{
  if (!fDim0 || fK11IsReduced)  {
    Simetrize();
    return (*fK11);
  }
  if (fK00->IsDecomposed()) {
    if( !fK01IsComputed ){
      Substitution((fK01));
      fK01IsComputed=1;
    }
  }
  else{
    Simetrize();
    fK00->SolveDirect((*fK01),fDecomposeType);
    fK01IsComputed=1;
  }


  //make [K11]=[k11]-[k10][k01]
  fK10->MultAdd((*fK01),(*fK11),(*fK11),-1,1);
  fK11IsReduced=1;
  return (*fK11);


}

void
TPZMatRed::U1(TPZFMatrix & F)
{

  K11Red();
  F1Red();
  F=(*fF1);
  fK11->SolveDirect( F ,ELU);


}

void
TPZMatRed::UGlobal(const TPZFMatrix & U1, TPZFMatrix & result)
{
  //[u0]=[A00^-1][F0]-[A00^-1][A01]
  if ( fK00->IsDecomposed() ){
    if( !fF0IsComputed ){
      //compute [F0]=[A00^-1][F0]
      Substitution((fF0));
      fF0IsComputed=1;
    }
    if( !fK01IsComputed ){
      //compute [A01]=[A00^-1][A01]
      Substitution( (fK01) );
      fK01IsComputed=1;
    }
  }else{
    Simetrize();
    fK00->SolveDirect( (*fF0),fDecomposeType );
    fF0IsComputed=1;
    Substitution( (fK01) );
    fK01IsComputed=1;
  }

  //make [u0]=[F0]-[U1]
  TPZFMatrix u0( fF0->Rows() , fF0->Cols() );
  fK01->MultAdd(U1,(*fF0),u0,-1,1);
  //compute result

  result.Redim( Rows(),fF0->Cols() );
  int c,r,r1;

  for(c=0; c<fF0->Cols(); c++){
    r1=0;
    for(r=0; r<fDim0; r++){
      result.PutVal( r,c,u0.GetVal(r,c) ) ;
    }
    //aqui r=fDim0
    for( ;r<Rows(); r++){
      result.PutVal( r,c,U1.GetVal(r1++,c) );
    }


  }


}

void
TPZMatRed::Print(const char *name , ostream &out ,const MatrixOutputFormat form ) const
{


  if(form != EInputFormat) {
    out << "Writing matrix 'TPZMatRed::" << name;
    out << "' (" << Dim() << " x " << Dim() << "):\n";

    fK00->Print("K(0,0)",out,form);
    fK01->Print("K(0,1)",out,form);
    fK10->Print("K(1,0)",out,form);
    fK11->Print("K(1,1)",out,form);


    if (fF0) {
      fF0->Print("F(0)",out,form);
      fF1->Print("F(1)",out,form);
    }


    out << "\n\n";
  } else {
    TPZMatrix::Print(name,out,form);
  }

}

int TPZMatRed::Substitution(TPZFMatrix *B) const{
	switch(fDecomposeType) {
		case ENoDecompose:
			default:
				Error( "TPZMatRed::Substitution called without initialized solver\n");
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

int TPZMatRed::Redim(int dim, int dim00){
	if(dim<dim00) Error("dim k00> dim");
	if(fK00) fK00->Redim(dim00,dim00);
	if(fK01) delete fK01;
	if(fK10) delete fK10;
	if(fK11) delete fK11;
	if(fF0)  delete fF0;
	if(fF1)  delete fF1;

	fDim0=dim00;
	fDim1=dim-dim00;
	fDecomposeType=ENoDecompose;
	fF0IsComputed=0;
	fK01IsComputed=0;
	fK11IsReduced=0;
	fF1IsReduced=0;

	fK01=new TPZFMatrix(fDim0,fDim1,0.);
	fK10=new TPZFMatrix(fDim1,fDim0,0.);
	fK11=new TPZFMatrix(fDim1,fDim1,0.);

	fF0=NULL;
	fF1=NULL;
	fRow = dim;
	fCol = dim;
	return 0;
}


int TPZMatRed::Zero(){
	if(fK00) fK00->Zero();
	if(fK01) fK01->Zero();
	if(fK10) fK10->Zero();
	if(fK11) fK11->Zero();
	if(fF0)  fF0->Zero();
	if(fF1)  fF1->Zero();
	fF0IsComputed=0;
	fK01IsComputed=0;
	fK11IsReduced=0;
	fF1IsReduced=0;
	return 0;
}
/*************************** Private ***************************/

/*************/
/*** Error ***/

int
TPZMatRed::Error(const char *msg ,const char *msg2) const
{
  cout << "TPZMatRed::" << msg << msg2 << ".\n";

  exit( 1 );
  return 0;
}
