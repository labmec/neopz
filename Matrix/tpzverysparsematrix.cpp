/**
 * @file
 * @brief Contains the implementation of the TPZVerySparseMatrix methods.
 */

#include "tpzverysparsematrix.h"

using namespace std;

template<class TVar>
TPZVerySparseMatrix<TVar>::TPZVerySparseMatrix() : fExtraSparseData()
{
}

template<class TVar>
TPZVerySparseMatrix<TVar>::~TPZVerySparseMatrix()
{
}

template<class TVar>
int TPZVerySparseMatrix<TVar>::PutVal(int row, int col, const TVar &val)
{
    if (row < 0 || col < 0 || row >this->fRow || col >this->fCol)
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

template<class TVar>
TPZVerySparseMatrix<TVar>::TPZVerySparseMatrix(const TPZFMatrix<TVar> &cp) : TPZMatrix<TVar>(cp)
{
	for(int i=0; i<this->fRow; i++)
	{
		for(int j=0; j<this->fCol; j++)
		{
			TVar a = cp.GetVal(i,j);
			if(a) PutVal(i,j,a);
		}
	}
}

template<class TVar>
const TVar & TPZVerySparseMatrix<TVar>::GetVal(int row, int col) const
{
    if (row < 0 || col < 0 || row >this->fRow || col >this->fCol)
    {
        cout<< "ERRO! em TPZVerySparseMatrix::GetVal: The row i or column j are incompatible with the rows or columns of a matrix"<< endl;
        this->gZero = 0.;
		return this->gZero;
    }
	
    pair <int,int> position(row,col);
    map<pair<int,int>, REAL>::const_iterator it;
    it = fExtraSparseData.find(position);
	
    if (it == fExtraSparseData.end() )
    {
        this->gZero = 0.;
        return this->gZero;
    }
    
    return it->second;
}

template<class TVar>
void TPZVerySparseMatrix<TVar>::MultAdd(const TPZFMatrix<TVar> & x, const TPZFMatrix<TVar> & y, TPZFMatrix<TVar> & z,
								  const TVar alpha, const TVar beta, const int opt, const int stride ) const
{
    if (!opt) 
    {
        if(this->Cols() != x.Rows()*stride || this->Rows() != y.Rows()*stride)
        {
            cout << "\nERROR! em TPZVerySparseMatrix::MultiplyAdd: incompatible dimensions in opt=false\n";
            return;
        } 
    } 
	else
		if (this->Rows() != x.Rows()*stride || this->Cols() != y.Rows()*stride)
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
    this->PrepareZ(y,z,beta,opt,stride);
    TVar val = 0.;
	
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
				TVar matrixval = it->second;
				
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
template<class TVar>
void TPZVerySparseMatrix<TVar>::Write(TPZStream &buf, int withclassid)
{
	TPZSaveable::Write(buf, withclassid);
	buf.Write(&this->fCol, 1);
	buf.Write(&this->fDecomposed, 1);
	buf.Write(&this->fDefPositive, 1);
	buf.Write(&this->fRow, 1);
	buf.Write(&this->gZero, 1);
	WriteMap(buf, withclassid, this->fExtraSparseData);
	
}
template<class TVar>
void TPZVerySparseMatrix<TVar>::WriteMap(TPZStream &buf, int withclassid, std::map<std::pair<int, int>, REAL> & TheMap)
{
	int mapsz = TheMap.size();
	buf.Write(&mapsz, 1);
	std::map<std::pair<int, int>, REAL>::iterator it;
	for(it = TheMap.begin(); it != TheMap.end(); it++)
	{
		int ii = 0, jj = 0;
		ii = it->first.first;
		jj = it->first.second;
		buf.Write(&ii, 1);
		buf.Write(&jj, 1);
		buf.Write(&it->second, 1);
	}
}

template<class TVar>
void TPZVerySparseMatrix<TVar>::Read(TPZStream &buf, void *context)
{
	TPZSaveable::Read(buf, context);
	buf.Read(&this->fCol, 1);
	buf.Read(&this->fDecomposed, 1);
	buf.Read(&this->fDefPositive, 1);
	buf.Read(&this->fRow, 1);
	buf.Read(&this->gZero, 1);
	ReadMap(buf, context, this->fExtraSparseData);
	
}
template<class TVar>
void TPZVerySparseMatrix<TVar>::ReadMap(TPZStream &buf, void *context, std::map<std::pair<int, int>, REAL> & TheMap)
{
	TheMap.clear();
	int size = 0;
	buf.Read(&size, 1);
	int i;
	for(i = 0; i < size; i++)
	{
		int ii = 0, jj = 0;
		TVar value = 0.;
		buf.Read(&ii, 1);
		buf.Read(&jj, 1);
		std::pair<int, int> item(ii, jj);
		buf.Read(&value, 1);
		std::pair<std::pair<int, int>, REAL > fullitem(item, value);
		TheMap.insert(fullitem);
	}
}

template class TPZVerySparseMatrix<REAL>;


