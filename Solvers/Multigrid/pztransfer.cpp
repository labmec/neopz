/**
 * @file
 * @brief Contains the implementation of the TPZTransfer methods. 
 */

#include "pztransfer.h"
#include "pzfmatrix.h"
#include <stdlib.h>

using namespace std;

template<class TVar>
TPZTransfer<TVar>::TPZTransfer() : TPZRegisterClassId(&TPZTransfer::ClassId),
fNTVarVar(0),fRowBlock(),fColBlock(),fColPosition(0),fNumberofColumnBlocks(0),fColumnBlockPosition(0),
fColumnBlockNumber(0),fColumnBlockLastUsed(0),fDoubleValues(0),fDoubleValLastUsed(0) {
}

// the sparse matrix blocks are defined by row, col
template<class TVar>
TPZTransfer<TVar>::TPZTransfer(TPZBlock &row, TPZBlock &col,int nvar, int nrowblocks, int ncolblocks) :
                                                        TPZRegisterClassId(&TPZTransfer::ClassId),
                                                        TPZMatrix<TVar>(), fNTVarVar(nvar), fRowBlock(), fColBlock(),
                                                        fColPosition(), fNumberofColumnBlocks(),
                                                        fColumnBlockPosition(0),fColumnBlockNumber(0),
                                                        fColumnBlockLastUsed(0),fDoubleValues(0),fDoubleValLastUsed(0)
{
	SetBlocks(row,col,nvar,nrowblocks,ncolblocks);
}

template<class TVar>
void TPZTransfer<TVar>::Print(const char *name,ostream &out,const MatrixOutputFormat form) const {
	if(form == EFormatted) {
		if(name) {
			out << name << endl;
		} else {
			out << "TPZTransfer<TVar>::Print : ";
		}
		out << "rows : " << this->Rows() << " cols : " << this->Cols() << endl;
	} else {
		out << this->Rows() << ' ' << this->Cols() << endl;
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
			TVar *locval = &fDoubleValues[fColumnBlockPosition[colpos+icbcounter]];
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
			TPZFMatrix<TVar> loc(rowsize,colsize,locval,rowsize*colsize);
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

//void TPZTransfer<TVar>::SetBlocks(TPZBlock &row,TPZBlock &col,int nvar, int nrowblocks, int ncolblocks){
template<class TVar>
void TPZTransfer<TVar>::SetBlocks(TPZBlock &row,TPZBlock &col,int nvar, int nrowblocks, int ncolblocks){
	// this operation will reset the matrix to zero
	// with no rows defined
	fRowBlock = row;
	fColBlock = col;
	fRowBlock.SetNBlocks(nrowblocks);
	fColBlock.SetNBlocks(ncolblocks);
	fNTVarVar = nvar;
	if(nvar!=1) {
		int i;
		for(i=0;i<nrowblocks;i++) fRowBlock.Set(i,row.Size(i)/nvar);
		for(i=0;i<ncolblocks;i++) fColBlock.Set(i,col.Size(i)/nvar);
		fRowBlock.Resequence();
		fColBlock.Resequence();
	}
	this->fRow = fRowBlock.Dim()*nvar;
	this->fCol = fColBlock.Dim()*nvar;
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

template<class TVar>
int TPZTransfer<TVar>::HasRowDefinition(int row){
	// returns 1 if the row is defined (i.e. has column entries)
	if(fColPosition[row] == -1) return 0;
	return 1;
}

template<class TVar>
void TPZTransfer<TVar>::AddBlockNumbers(int row, TPZVec<int> &colnumbers){
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

template<class TVar>
void TPZTransfer<TVar>::ExpandColumnVectorEntries(int num){
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


template<class TVar>
void TPZTransfer<TVar>::SetBlockMatrix(int row, int col, TPZFMatrix<TVar> &mat){
	// sets the row,col block equal to matrix mat
	// if row col was not specified by AddBlockNumbers, an error
	//		will be issued and exit
	// find row,col
	int colpos = fColPosition[row];
	int numcolblocks = fNumberofColumnBlocks[row];
	if(colpos == -1 || numcolblocks == -1) {
		cout << "TPZTransfer<TVar>::SetBlockMatrix called for ilegal parameters : "
		<< " row = " << row << " col = " << col << " colpos = " << colpos <<
		" numcolblocks = " << numcolblocks << endl;
		DebugStop();
	}
	int ic = colpos, lastic = colpos+numcolblocks;
	for(;ic<lastic;ic++) {
		if(fColumnBlockNumber[ic] == col) break;
	}
	if(ic == lastic) {
		cout << "TPZTransfer<TVar>::SetBlockMatrix column not found for row = " << row <<
		" col = " << col << endl;
		DebugStop();
	}
	int nblrows = fRowBlock.Size(row);
	int nblcols = fColBlock.Size(col);
	if(nblrows != mat.Rows() || nblcols != mat.Cols()) {
		cout << "TPZTransfer<TVar>::SetBlockMatrix matrix has incompatible dimensions : "
		" nblrows = " << nblrows << " nblcols = " << nblcols << " mat.rows = " <<
		mat.Rows() << " mat.cols " << mat.Cols() << endl;
		DebugStop();
	}
	ExpandDoubleValueEntries(nblrows*nblcols);
	TPZFMatrix<TVar> bl(nblrows,nblcols,&fDoubleValues[fDoubleValLastUsed],nblrows*nblcols);
	bl = mat;
	fColumnBlockPosition[ic] = fDoubleValLastUsed;
	fDoubleValLastUsed += nblrows*nblcols;
}

template<class TVar>
void TPZTransfer<TVar>::ExpandDoubleValueEntries(int num){
	if(fDoubleValLastUsed+num > fDoubleValues.NAlloc()) {
		fDoubleValues.Expand(fDoubleValues.NAlloc()+4000);
		//    Print("After the double value expansion");
	}
	int nextsize = fDoubleValues.NElements()+num;
	fDoubleValues.Resize(nextsize);
}


template<class TVar>
void TPZTransfer<TVar>::MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
                                      TVar alpha, TVar beta, int opt) const{
    // multiplies the transfer matrix and puts the result in z
    if ((!opt && this->Cols() != x.Rows()) || (opt && this->Rows() != x.Rows()))
        this->Error( "TPZTransfer<TVar>::MultAdd <matrices with incompatible dimensions>" );
    if(x.Cols() != y.Cols() || x.Cols() != z.Cols()) {
        this->Error ("TPZTransfer<TVar>::MultiplyAdd incompatible dimensions\n");
    }
    int rows = fRowBlock.MaxBlockSize();
    int xcols = x.Cols();
    int ic, c, r;
    this->PrepareZ(y,z,beta,opt);
    if (fNTVarVar == 1) {
        MultAddScalar(x, y, z, alpha, beta, opt);
    }
    else
    {
        int nrc = x.Rows();
        int ncc = x.Cols();
        int thisr = this->Rows()/fNTVarVar;
        int thisc = this->Cols()/fNTVarVar;
        for(int iv=0; iv<fNTVarVar; iv++) {
            TPZFMatrix<TVar> tempcoarse(thisc,ncc), tempfine(thisr,ncc);
            for (int i=0; i<thisr; i++) {
                for (int c=0; c<ncc; c++) {
                    tempcoarse(i,c) = x.GetVal(iv+i*fNTVarVar,c);
                }
            }
            MultAddScalar(tempcoarse, tempfine, tempfine, alpha, 0., opt);
            for (int i=0; i<thisr; i++) {
                for (int c=0; c<ncc; c++) {
                    z(iv+i*fNTVarVar,c) = tempfine.GetVal(i,c);
                }
            }
        }
    }
}

template<class TVar>
void TPZTransfer<TVar>::MultAddScalar(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
						  TVar alpha, TVar beta, int opt) const{
	// multiplies the transfer matrix and puts the result in z
	if ((!opt && this->Cols() != x.Rows()) || (opt && this->Rows() != x.Rows()))
		this->Error( "TPZTransfer<TVar>::MultAdd <matrices with incompatible dimensions>" );
	if(x.Cols() != y.Cols() || x.Cols() != z.Cols()) {
		this->Error ("TPZTransfer<TVar>::MultiplyAdd incompatible dimensions\n");
	}
	int rows = fRowBlock.MaxBlockSize();
	int xcols = x.Cols();
	int ic, c, r;
	this->PrepareZ(y,z,beta,opt);
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
				TPZFMatrix<TVar> zloc(rowblocksize,1,&z(rowblockpos,ic),rowblocksize);
				for( c=0; c < numcolbl; c++) {
					int col = fColumnBlockNumber[colpos+c];
					int colblockpos = fColBlock.Position(col);
					int colblocksize = fColBlock.Size(col);
					TPZFMatrix<TVar> xloc(colblocksize,1,&x.g(colblockpos,ic), colblocksize);
					TPZFMatrix<TVar> aloc(rowblocksize,colblocksize,
									&fDoubleValues[fColumnBlockPosition[colpos+c]],rowblocksize*colblocksize);
					aloc.MultAdd(xloc,zloc,zloc,alpha,1.,opt);
				}
			}
		} else {
			for ( r=0; r < rows; r++) {
				int rowblockpos = fRowBlock.Position(r);
				int rowblocksize = fRowBlock.Size(r);
				if(!rowblocksize) continue;
				TPZFMatrix<TVar> xloc(rowblocksize,1,&x.g(rowblockpos,ic),rowblocksize);
				int colpos = fColPosition[r];
				if(colpos == -1) continue;
				int numcolbl = fNumberofColumnBlocks[r];
				if(numcolbl == 0) continue;
				for( c=0; c < numcolbl; c++) {
					int col = fColumnBlockNumber[colpos+c];
					int colblockpos = fColBlock.Position(col);
					int colblocksize = fColBlock.Size(col);
					TPZFMatrix<TVar> zloc(colblocksize,1,&z(colblockpos,ic),colblocksize);
					TPZFMatrix<TVar> aloc(rowblocksize,colblocksize,
									&fDoubleValues[fColumnBlockPosition[colpos+c]],rowblocksize*colblocksize);
					aloc.MultAdd(xloc,zloc,zloc,alpha,1.,opt);
					col = rowblockpos-1;
					rowblockpos = col+1;
				}
			}
		}
	}
}

/**
 * Will transfer the solution, taking into acount there may be more than
 * one TVar variable
 */
template<class TVar>
void TPZTransfer<TVar>::TransferSolution(const TPZFMatrix<TVar> &coarsesol, TPZFMatrix<TVar> &finesol){
	int iv;
	int nrf = finesol.Rows();
	int ncf = finesol.Cols();
	int nrc = coarsesol.Rows();
	int ncc = coarsesol.Cols();
    int thisr = this->Rows();
    int thisc = this->Cols();
	if(nrf != this->Rows()*fNTVarVar || ncf != ncc) {
		nrf = this->Rows()*fNTVarVar;
		ncf = ncc;
		finesol.Redim(nrf,ncf);
	}
    if (fNTVarVar == 1) {
        MultiplyScalar(coarsesol,finesol,0);
    }
    else
    {
        for(iv=0; iv<fNTVarVar; iv++) {
            TPZFMatrix<TVar> tempcoarse(thisc,ncc), tempfine(thisr,ncf);
            for (int i=0; i<thisr; i++) {
                for (int c=0; c<ncf; c++) {
                    tempcoarse(i,c) = coarsesol.GetVal(iv+i*fNTVarVar,c);

                }
            }
            MultiplyScalar(tempcoarse, tempfine, 0);
            for (int i=0; i<thisr; i++) {
                for (int c=0; c<ncf; c++) {
                    finesol(iv+i*fNTVarVar,c) = tempfine.GetVal(i,c);
                    
                }
            }
        }
    }
}

/**
 * Will transfer the residual, taking into acount there may be more than
 * one TVar variable
 */
template<class TVar>
void TPZTransfer<TVar>::TransferResidual(const TPZFMatrix<TVar> &fine, TPZFMatrix<TVar> &coarse){
	int iv;
	int nrf = fine.Rows();
	int ncf = fine.Cols();
	int nrc = coarse.Rows();
	int ncc = coarse.Cols();
    int thisr = this->Rows();
    int thisc = this->Cols();
	if(ncc != ncf || nrc != this->Cols()*fNTVarVar) {
		ncc = ncf;
		nrc = this->Cols()*fNTVarVar;
		coarse.Redim(nrc,ncf);
	}
    if (fNTVarVar == 1) {
        MultiplyScalar(fine,coarse,1);
    }
    else
    {
        for(iv=0; iv<fNTVarVar; iv++) {
            TPZFMatrix<TVar> tempcoarse(thisr,ncc), tempfine(thisc,ncf);
            for (int i=0; i<thisc; i++) {
                for (int c=0; c<ncf; c++) {
                    tempfine(i,c) = fine.GetVal(iv+i*fNTVarVar,c);
                    
                }
            }
            MultiplyScalar(tempfine, tempcoarse, 1);
            for (int i=0; i<thisr; i++) {
                for (int c=0; c<ncf; c++) {
                    coarse(iv+i*fNTVarVar,c) = tempcoarse(i,c);
                    
                }
            }
        }
    }

}

template<class TVar>
void TPZTransfer<TVar>::MultiplyScalar(const TPZFMatrix<TVar> &A, TPZFMatrix<TVar> &B,const int opt) const {
	if ((opt==0 && this->Cols() != A.Rows()) || (opt ==1 && this->Rows() != A.Rows()))
    {
		this->Error( "TPZTransfer<TVar>::Multiply incompatible dimensions" );
	}
	if(!opt && (B.Rows()*fNTVarVar != this->Rows() || B.Cols() != A.Cols())) {
		B.Redim(this->Rows()/fNTVarVar,A.Cols()/fNTVarVar);
	}
	else if (opt && (B.Rows()*fNTVarVar != this->Cols() || B.Cols() != A.Cols())) {
		B.Redim(this->Cols()/fNTVarVar,A.Cols()/fNTVarVar);
	}
	MultAddScalar( A, B, B, 1.0, 0.0, opt);
}

template class TPZTransfer<float>;
template class TPZTransfer<double>;
template class TPZTransfer<long double>;

template class TPZTransfer<std::complex<double> >;
template class TPZTransfer<std::complex<float> >;
template class TPZTransfer<std::complex<long double> >;

