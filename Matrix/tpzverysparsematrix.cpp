
#include "tpzverysparsematrix.h"

using namespace std;

TPZVerySparseMatrix::TPZVerySparseMatrix() : fExtraSparseData()
{
}

TPZVerySparseMatrix::~TPZVerySparseMatrix()
{
}

int TPZVerySparseMatrix::PutVal(int row, int col, const REAL &val)
{
    if (row < 0 || col < 0 || row >fRow || col >fCol)
    {
        cout<< "ERRO! em TPZVerySparseMatrix::PutVal: The row i or column j are incompatible with the rows or columns of a matrix"<< endl;
        return -1;
   }
  
    pair <int,int> position(row,col);
	std::map <std::pair<int, int>, REAL>::iterator it = fExtraSparseData.find(position);
	if(val == 0. && it != fExtraSparseData.end()) 
	{
		fExtraSparseData.erase(it);
	}
	else if(val && it != fExtraSparseData.end())
	{
		it->second = val;
	}
	else if (val)
	{
		fExtraSparseData[position] = val;
	}
	return 0;
}

TPZVerySparseMatrix::TPZVerySparseMatrix(const TPZFMatrix &cp) : TPZMatrix(cp)
{
	for(int i=0; i<fRow; i++)
	{
		for(int j=0; j<fCol; j++)
		{
			REAL a = cp.GetVal(i,j);
			if(a) PutVal(i,j,a);
		}
	}
}


const REAL & TPZVerySparseMatrix::GetVal(int row, int col) const
{
    if (row < 0 || col < 0 || row >fRow || col >fCol)
    {
        cout<< "ERRO! em TPZVerySparseMatrix::GetVal: The row i or column j are incompatible with the rows or columns of a matrix"<< endl;
        gZero = 0.;
		return gZero;
    }
  
    pair <int,int> position(row,col);
    map<pair<int,int>, REAL>::const_iterator it;
    it = fExtraSparseData.find(position);
  
    if (it == fExtraSparseData.end() )
    {
        gZero = 0.;
        return gZero;
    }
    
    return it->second;
}

void TPZVerySparseMatrix::MultAdd(const TPZFMatrix & x, const TPZFMatrix & y, TPZFMatrix & z,
                                                        const REAL alpha, const REAL beta, const int opt, const int stride ) const
{
    if (!opt) 
    {
        if(Cols() != x.Rows()*stride || Rows() != y.Rows()*stride)
        {
            cout << "\nERROR! em TPZVerySparseMatrix::MultiplyAdd: incompatible dimensions in opt=false\n";
            return;
        } 
    } 
      else
          if (Rows() != x.Rows()*stride || Cols() != y.Rows()*stride)
          {
             cout << "\nERROR! em TPZVerySparseMatrix::MultiplyAdd: incompatible dimensions in opt=true\n";
              return; 
          }
    if(x.Cols() != y.Cols() || x.Cols() != z.Cols() || y.Rows() != z.Rows() )
    {
        cout << "\nERROR! em TPZVerySparseMatrix::MultiplyAdd : incompatible dimensions in x, y or z\n";
        return;
    }
    
    int xcols = x.Cols();
    int ic, c, r;
    PrepareZ(y,z,beta,opt,stride);
    REAL val = 0.;
  
    for (ic = 0; ic < xcols; ic++)
    {
        if(!opt) 
        {
          map< pair<int,int>, REAL>::const_iterator it;
          
          for(it = fExtraSparseData.begin(); it!= fExtraSparseData.end(); it++)
          {
              pair <int, int> position(it->first);
              c = position.second;
              r = position.first;
              REAL matrixval = it->second;
      
              val = z(r*stride,ic) + alpha * matrixval * x.GetVal(c*stride,ic);
              z.PutVal(r*stride,ic,val);
          }
        } else 
          {
              map<pair<int,int>, REAL>::const_iterator it;
              
              for(it = fExtraSparseData.begin(); it != fExtraSparseData.end(); it++)
              {
                  pair <int, int> posicao(it->first);
                  r = posicao.second;
                  c = posicao.first;
                  double matrixval = it->second;
                  z(r*stride,ic) += (alpha*matrixval)*x.GetVal(c*stride,ic);
              }
          }
    }
}

