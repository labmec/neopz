/**
 * @file
 * @brief Contains the implementation of the TPZTransfer methods. 
 */
#include "pztransfer.h"
#include "pzfmatrix.h"
#include <stdlib.h>

using namespace std;

TPZTransfer::TPZTransfer() :
fNStateVar(0),fRowBlock(),fColBlock(),fColPosition(0),fNumberofColumnBlocks(0),fColumnBlockPosition(0),
fColumnBlockNumber(0),fColumnBlockLastUsed(0),fDoubleValues(0),fDoubleValLastUsed(0) {
}

TPZTransfer::TPZTransfer(TPZBlock<REAL> &row, TPZBlock<REAL> &col,int nvar, int nrowblocks, int ncolblocks) :
// the sparse matrix blocks are defined by row, col
TPZMatrix<REAL>(), fNStateVar(nvar), fRowBlock(), fColBlock(),
fColPosition(), fNumberofColumnBlocks(),
fColumnBlockPosition(0),fColumnBlockNumber(0),
fColumnBlockLastUsed(0),fDoubleValues(0),fDoubleValLastUsed(0) {
	SetBlocks(row,col,nvar,nrowblocks,ncolblocks);
}

void TPZTransfer::Print(const char *name,ostream &out,const MatrixOutputFormat form) const {
	if(form == EFormatted) {
		if(name) {
			out << name << endl;
		} else {
			out << "TPZTransfer::Print : ";
		}
		out << "rows : " << Rows() << " cols : " << Cols() << endl;
	} else {
		out << Rows() << ' ' << Cols() << endl;
	}
	int irb;
	for(irb=0 ; irb < fRowBlock.MaxBlockSize(); irb++) {
		if(form == EFormatted) {
			out << "block row number : " << irb << endl;
			out << "number of column blocks : " << fNumberofColumnBlocks[irb] << endl;
			out << "position of first column block : " << fColPosition[irb] << endl;
		}
		int colpos = fColPosition[irb];
		int numcolbl = fNumberofColumnBlocks[irb];
		int icbcounter;
		for(icbcounter=0; icbcounter < numcolbl; icbcounter++) {
			int icb = fColumnBlockNumber[colpos+icbcounter];
			REAL *locval = &fDoubleValues[fColumnBlockPosition[colpos+icbcounter]];
			if(form == EFormatted) {
				out << "column block counter : " << icbcounter <<
				" column block number " << icb << endl;
			}
			int rownumber = fRowBlock.Position(irb);
			int colnumber = fColBlock.Position(icb);
			int rowsize = fRowBlock.Size(irb);
			int colsize = fColBlock.Size(icb);
			if(form == EFormatted) {
				out << "row position number : " << rownumber <<
				" column position number " << colnumber << endl;
				out << "block sizes : row : " << rowsize << " col "  << colsize << endl;
			}
			TPZFMatrix<REAL> loc(rowsize,colsize,locval,rowsize*colsize);
			if(form == EFormatted) {
				loc.Print(0,out,form);
			} else {
				int ir,ic;
				for(ir=0; ir<rowsize; ir++) {
					for(ic=0; ic<colsize; ic++) {
						out << (rownumber+ir) << ' ' << (colnumber+ic) << ' ' << loc(ir,ic) << endl;
					}
				}
			}
		}
	}
}

void TPZTransfer::SetBlocks(TPZBlock<REAL> &row,TPZBlock<REAL> &col,int nvar, int nrowblocks, int ncolblocks){
	// this operation will reset the matrix to zero
	// with no rows defined
	fRowBlock = row;
	fColBlock = col;
	fRowBlock.SetNBlocks(nrowblocks);
	fColBlock.SetNBlocks(ncolblocks);
	fNStateVar = nvar;
	if(nvar!=1) {
		int i;
		for(i=0;i<nrowblocks;i++) fRowBlock.Set(i,row.Size(i)/nvar);
		for(i=0;i<ncolblocks;i++) fColBlock.Set(i,col.Size(i)/nvar);
		fRowBlock.Resequence();
		fColBlock.Resequence();
	}
	fRow = fRowBlock.Dim();
	fCol = fColBlock.Dim();
	fColPosition.Resize(row.MaxBlockSize());
	fColPosition.Fill(-1,0);
	fNumberofColumnBlocks.Resize(row.MaxBlockSize());
	fNumberofColumnBlocks.Fill(-1,0);
	fColumnBlockPosition.Resize(0);
	fColumnBlockNumber.Resize(0);
	fColumnBlockLastUsed = 0;
	fDoubleValues.Resize(0);
	fDoubleValLastUsed = 0;
}

int TPZTransfer::HasRowDefinition(int row){
	// returns 1 if the row is defined (i.e. has column entries)
	if(fColPosition[row] == -1) return 0;
	return 1;
}

void TPZTransfer::AddBlockNumbers(int row, TPZVec<int> &colnumbers){
	if(HasRowDefinition(row)) {
		cout << "TPZTransfer:SetBlocks called for an already defined row = " <<
		row << endl;
		DebugStop();
	}
	fColPosition[row] = fColumnBlockLastUsed;
	fNumberofColumnBlocks[row] = colnumbers.NElements();
	// will specify the sparsity pattern of row
	ExpandColumnVectorEntries(colnumbers.NElements());
	int ic= fColumnBlockLastUsed, lastic = ic+colnumbers.NElements();
	for(;ic < lastic; ic++) {
		fColumnBlockNumber[ic] = colnumbers[ic-fColumnBlockLastUsed];
	}
	fColumnBlockLastUsed = lastic;
}

void TPZTransfer::ExpandColumnVectorEntries(int num){
	while(fColumnBlockLastUsed+num > fColumnBlockNumber.NAlloc()) {
		int nextsize = fColumnBlockNumber.NAlloc()+1000;
		fColumnBlockNumber.Expand(nextsize);
		fColumnBlockPosition.Expand(nextsize);
	}
	fColumnBlockNumber.Resize(fColumnBlockLastUsed+num);
	fColumnBlockNumber.Fill(-1,fColumnBlockLastUsed);
	fColumnBlockPosition.Resize(fColumnBlockLastUsed+num);
	fColumnBlockPosition.Fill(0,fColumnBlockLastUsed);
}


void TPZTransfer::SetBlockMatrix(int row, int col, TPZFMatrix<REAL> &mat){
	// sets the row,col block equal to matrix mat
	// if row col was not specified by AddBlockNumbers, an error
	//		will be issued and exit
	// find row,col
	int colpos = fColPosition[row];
	int numcolblocks = fNumberofColumnBlocks[row];
	if(colpos == -1 || numcolblocks == -1) {
		cout << "TPZTransfer::SetBlockMatrix called for ilegal parameters : "
		<< " row = " << row << " col = " << col << " colpos = " << colpos <<
		" numcolblocks = " << numcolblocks << endl;
		DebugStop();
	}
	int ic = colpos, lastic = colpos+numcolblocks;
	for(;ic<lastic;ic++) {
		if(fColumnBlockNumber[ic] == col) break;
	}
	if(ic == lastic) {
		cout << "TPZTransfer::SetBlockMatrix column not found for row = " << row <<
		" col = " << col << endl;
		DebugStop();
	}
	int nblrows = fRowBlock.Size(row);
	int nblcols = fColBlock.Size(col);
	if(nblrows != mat.Rows() || nblcols != mat.Cols()) {
		cout << "TPZTransfer::SetBlockMatrix matrix has incompatible dimensions : "
		" nblrows = " << nblrows << " nblcols = " << nblcols << " mat.rows = " <<
		mat.Rows() << " mat.cols " << mat.Cols() << endl;
		DebugStop();
	}
	ExpandDoubleValueEntries(nblrows*nblcols);
	TPZFMatrix<REAL> bl(nblrows,nblcols,&fDoubleValues[fDoubleValLastUsed],nblrows*nblcols);
	bl = mat;
	fColumnBlockPosition[ic] = fDoubleValLastUsed;
	fDoubleValLastUsed += nblrows*nblcols;
}

void TPZTransfer::ExpandDoubleValueEntries(int num){
	if(fDoubleValLastUsed+num > fDoubleValues.NAlloc()) {
		fDoubleValues.Expand(fDoubleValues.NAlloc()+4000);
		//    Print("After the double value expansion");
	}
	int nextsize = fDoubleValues.NElements()+num;
	fDoubleValues.Resize(nextsize);
}


void TPZTransfer::MultAdd(const TPZFMatrix<REAL> &x,const TPZFMatrix<REAL> &y, TPZFMatrix<REAL> &z,
						  REAL alpha, REAL beta, int opt, int stride) const{
	// multiplies the transfer matrix and puts the result in z
	if ((!opt && Cols()*stride != x.Rows()) || (opt && Rows()*stride != x.Rows()))
		Error( "TPZTransfer::MultAdd <matrices with incompatible dimensions>" );
	if(x.Cols() != y.Cols() || x.Cols() != z.Cols()) {
		Error ("TPZTransfer::MultiplyAdd incompatible dimensions\n");
	}
	int rows = fRowBlock.MaxBlockSize();
	int xcols = x.Cols();
	int ic, c, r;
	PrepareZ(y,z,beta,opt,stride);
	for (ic = 0; ic < xcols; ic++) {
		if(!opt) {
			for ( r=0; r < rows; r++) {
				int rowblockpos = fRowBlock.Position(r);
				int rowblocksize = fRowBlock.Size(r);
				if(!rowblocksize) continue;
				int colpos = fColPosition[r];
				if(colpos == -1) continue;
				int numcolbl = fNumberofColumnBlocks[r];
				if(numcolbl == 0) continue;
				TPZFMatrix<REAL> zloc(rowblocksize*stride,1,&z(rowblockpos*stride,ic),rowblocksize*stride);
				for( c=0; c < numcolbl; c++) {
					int col = fColumnBlockNumber[colpos+c];
					int colblockpos = fColBlock.Position(col);
					int colblocksize = fColBlock.Size(col);
					TPZFMatrix<REAL> xloc(colblocksize*stride,1,&x.g(colblockpos*stride,ic), colblocksize*stride);
					TPZFMatrix<REAL> aloc(rowblocksize,colblocksize,
									&fDoubleValues[fColumnBlockPosition[colpos+c]],rowblocksize*colblocksize);
					aloc.MultAdd(xloc,zloc,zloc,alpha,1.,opt,stride);
				}
			}
		} else {
			for ( r=0; r < rows; r++) {
				int rowblockpos = fRowBlock.Position(r);
				int rowblocksize = fRowBlock.Size(r);
				if(!rowblocksize) continue;
				TPZFMatrix<REAL> xloc(rowblocksize*stride,1,&x.g(rowblockpos*stride,ic),rowblocksize*stride);
				int colpos = fColPosition[r];
				if(colpos == -1) continue;
				int numcolbl = fNumberofColumnBlocks[r];
				if(numcolbl == 0) continue;
				for( c=0; c < numcolbl; c++) {
					int col = fColumnBlockNumber[colpos+c];
					int colblockpos = fColBlock.Position(col);
					int colblocksize = fColBlock.Size(col);
					TPZFMatrix<REAL> zloc(colblocksize*stride,1,&z(colblockpos*stride,ic),colblocksize*stride);
					TPZFMatrix<REAL> aloc(rowblocksize,colblocksize,
									&fDoubleValues[fColumnBlockPosition[colpos+c]],rowblocksize*colblocksize);
					aloc.MultAdd(xloc,zloc,zloc,alpha,1.,opt,stride);
					col = rowblockpos-1;
					rowblockpos = col+1;
				}
			}
		}
	}
}

/**
 * Will transfer the solution, taking into acount there may be more than
 * one state variable
 */
void TPZTransfer::TransferSolution(const TPZFMatrix<REAL> &coarsesol, TPZFMatrix<REAL> &finesol){
	int iv;
	int nrf = finesol.Rows();
	int ncf = finesol.Cols();
	int nrc = coarsesol.Rows();
	int ncc = coarsesol.Cols();
	if(nrf != Rows()*fNStateVar || ncf != ncc) {
		nrf = Rows()*fNStateVar;
		ncf = ncc;
		finesol.Redim(nrf,ncf);
	}
	for(iv=0; iv<fNStateVar; iv++) {
		REAL *fvp = &finesol.g(iv,0);
		REAL *cvp = &coarsesol.g(iv,0);
		TPZFMatrix<REAL> finewrap(nrf,ncf,fvp,nrf*ncf);
		TPZFMatrix<REAL> coarsewrap(nrc,ncc,cvp,nrc*ncc);
		Multiply(coarsewrap,finewrap,0,fNStateVar);
	}
}

/**
 * Will transfer the residual, taking into acount there may be more than
 * one state variable
 */

void TPZTransfer::TransferResidual(const TPZFMatrix<REAL> &fine, TPZFMatrix<REAL> &coarse){
	int iv;
	int nrf = fine.Rows();
	int ncf = fine.Cols();
	int nrc = coarse.Rows();
	int ncc = coarse.Cols();
	if(ncc != ncf || nrc != Cols()*fNStateVar) {
		ncc = ncf;
		nrc = Cols()*fNStateVar;
		coarse.Redim(nrc,ncf);
	}
	for(iv=0; iv<fNStateVar; iv++) {
		REAL *fvp = &fine.g(iv,0);
		REAL *cvp = &coarse.g(iv,0);
		TPZFMatrix<REAL> finewrap(nrf,ncf,fvp,nrf*ncf);
		TPZFMatrix<REAL> coarsewrap(nrc,ncc,cvp,nrc*ncc);
		Multiply(finewrap, coarsewrap,1,fNStateVar);
	}
}

void TPZTransfer::Multiply(const TPZFMatrix<REAL> &A, TPZFMatrix<REAL> &B,const int opt,
						   const int stride) const {
	if ((opt==0 && Cols()*stride != A.Rows()) || (opt ==1 && Rows()*stride != A.Rows()))
    {
		Error( "TPZTransfer::Multiply incompatible dimensions" );
	}
	if(!opt && (B.Rows() != Rows()*stride || B.Cols() != A.Cols())) {
		B.Redim(Rows()*stride,A.Cols());
	}
	else if (opt && (B.Rows() != Cols()*stride || B.Cols() != A.Cols())) {
		B.Redim(Cols()*stride,A.Cols());
	}
	//   if(opt == 0) {
	//     B.Redim(Rows()*stride, A.Cols() );
	//   } else {
	//     B.Redim(Cols()*stride, A.Cols() );
	//   }
	
	MultAdd( A, B, B, 1.0, 0.0, opt,stride);
}
