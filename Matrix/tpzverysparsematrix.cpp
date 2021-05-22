/**
 * @file
 * @brief Contains the implementation of the TPZVerySparseMatrix methods.
 */

#include "tpzverysparsematrix.h"
#include "pzfmatrix.h"
#include "TPZStream.h"
using namespace std;

template<class TVar>
TPZVerySparseMatrix<TVar>::TPZVerySparseMatrix() : 
TPZRegisterClassId(&TPZVerySparseMatrix::ClassId),
fExtraSparseData()
{
}

template<class TVar>
TPZVerySparseMatrix<TVar>::~TPZVerySparseMatrix()
{
}

template<class TVar>
int TPZVerySparseMatrix<TVar>::PutVal(const int64_t row,const int64_t col, const TVar &val)
{
    if (row < 0 || col < 0 || row >this->fRow || col >this->fCol)
    {
        cout<< "ERRO! em TPZVerySparseMatrix::PutVal: The row i or column j are incompatible with the rows or columns of a matrix"<< endl;
        return -1;
	}
	
    pair <int64_t,int64_t> position(row,col);
	typename std::map <std::pair<int64_t, int64_t>, TVar>::iterator it = this->fExtraSparseData.find(position);
	if(val == (TVar)0. && it != fExtraSparseData.end()) 
	{
		fExtraSparseData.erase(it);
	}
	else if( it != fExtraSparseData.end())
	{
		it->second = val;
	}
	else
	{
		fExtraSparseData[position] = val;
	}
	return 0;
}

template<class TVar>
TPZVerySparseMatrix<TVar>::TPZVerySparseMatrix(const TPZFMatrix<TVar> &cp) : TPZRegisterClassId(&TPZVerySparseMatrix::ClassId),TPZMatrix<TVar>(cp)
{
	for(int64_t i=0; i<this->fRow; i++)
	{
		for(int64_t j=0; j<this->fCol; j++)
		{
			TVar a = cp.GetVal(i,j);
			if(!IsZero(a)) PutVal(i,j,a);
		}
	}
}

template<class TVar>
void TPZVerySparseMatrix<TVar>::Simetrize() {
  
  int64_t rows = this->Rows();
  int64_t cols = this->Cols();
  
  // TODO: introduce error handling mechanism  
  if ( rows != cols) {
      cerr << "Error: Simetrize only work for square matrices";
      return;
  }
  
  typename std::map <std::pair<int64_t, int64_t>, TVar>::iterator it = this->fExtraSparseData.begin();
  typename std::map <std::pair<int64_t, int64_t>, TVar>::iterator end = this->fExtraSparseData.end();
  typename std::map <std::pair<int64_t, int64_t>, TVar>::iterator next;
  
  std::list< std::pair< std::pair<int64_t, int64_t>, TVar > > temp;
  
  for(; it != end; it=next) {
    const std::pair<int64_t, int64_t>& key = it->first;
    next = it;
    next++;
    
    if(key.first < key.second) {
      temp.push_back( std::pair< std::pair<int64_t, int64_t>, TVar > (std::pair<int64_t, int64_t>(key.second, key.first), it->second));
    }
    else if (key.first > key.second) {
      this->fExtraSparseData.erase(it);
    }
  }
  
  typename std::list< std::pair< std::pair<int64_t, int64_t>, TVar > >::iterator at = temp.begin();
  typename std::list< std::pair< std::pair<int64_t, int64_t>, TVar > >::iterator atEnd = temp.end();
  
  	for(; at != atEnd; at++) {
		this->fExtraSparseData [ at->first ] = at->second;
	}

}

template<class TVar>
void TPZVerySparseMatrix<TVar>::Transpose(TPZVerySparseMatrix<TVar> *T) const {
  int64_t rows = this->Rows();
  int64_t cols = this->Cols();
  T->Resize( cols, rows );
  T->fExtraSparseData.clear();
  typename std::map <std::pair<int64_t, int64_t>, TVar>::const_iterator it = this->fExtraSparseData.begin();
  typename std::map <std::pair<int64_t, int64_t>, TVar>::const_iterator end = this->fExtraSparseData.end();

  for (; it != end; it++) {
    const std::pair<int64_t, int64_t>& key = it->first;
    T->fExtraSparseData[std::pair<int64_t,int64_t>(key.second, key.first)] = it->second;
  }
}

template<class TVar>
const TVar TPZVerySparseMatrix<TVar>::GetVal(const int64_t row, const int64_t col) const
{
    if (row < 0 || col < 0 || row >this->fRow || col >this->fCol)
    {
        cout<< "ERROR! em TPZVerySparseMatrix::GetVal: The row i or column j are incompatible with the rows or columns of a matrix"<< endl;
		return (TVar)0;
    }
	
    pair <int64_t,int64_t> position(row,col);
    typename map<pair<int64_t,int64_t>, TVar>::const_iterator it;
    it = fExtraSparseData.find(position);
	
    if (it == fExtraSparseData.end() )
    {
        return (TVar) 0;
    }
    
    return it->second;
}

template<class TVar>
TVar* &TPZVerySparseMatrix<TVar>::Elem()
{
  PZError<<__PRETTY_FUNCTION__;
  PZError<<"ERROR:Not implemented\n.Aborting...\n";
  DebugStop();
  static TVar* t{nullptr};
  return t;
}

template<class TVar>
const TVar* TPZVerySparseMatrix<TVar>::Elem()const
{
  PZError<<__PRETTY_FUNCTION__;
  PZError<<"ERROR:Not implemented\n.Aborting...\n";
  DebugStop();
  return nullptr;
}
template<class TVar>
void TPZVerySparseMatrix<TVar>::MultAdd(const TPZFMatrix<TVar> & x, const TPZFMatrix<TVar> & y, TPZFMatrix<TVar> & z,
								  const TVar alpha, const TVar beta, const int opt) const
{
    if (!opt) 
    {
        if(this->Cols() != x.Rows() || this->Rows() != y.Rows())
        {
            cout << "\nERROR! em TPZVerySparseMatrix::MultiplyAdd: incompatible dimensions in opt=false\n";
            return;
        } 
    } 
	else
		if (this->Rows() != x.Rows() || this->Cols() != y.Rows())
		{
			cout << "\nERROR! em TPZVerySparseMatrix::MultiplyAdd: incompatible dimensions in opt=true\n";
			return; 
		}
    if(x.Cols() != y.Cols() || x.Cols() != z.Cols() || y.Rows() != z.Rows() )
    {
        cout << "\nERROR! em TPZVerySparseMatrix::MultiplyAdd : incompatible dimensions in x, y or z\n";
        return;
    }
    
    int64_t xcols = x.Cols();
    int64_t ic, c, r;
    this->PrepareZ(y,z,beta,opt);
    TVar val = 0.;
	
    for (ic = 0; ic < xcols; ic++)
    {
        if(!opt) 
        {
			typename map< pair<int64_t,int64_t>, TVar>::const_iterator it;
			
			for(it = fExtraSparseData.begin(); it!= fExtraSparseData.end(); it++)
			{
				pair <int64_t, int64_t> position(it->first);
				c = position.second;
				r = position.first;
				TVar matrixval = it->second;
				
				val = z(r,ic) + alpha * matrixval * x.GetVal(c,ic);
				z.PutVal(r,ic,val);
			}
        } else 
		{
			typename map<pair<int64_t,int64_t>, TVar>::const_iterator it;
			
			for(it = fExtraSparseData.begin(); it != fExtraSparseData.end(); it++)
			{
				pair <int64_t, int64_t> posicao(it->first);
				r = posicao.second;
				c = posicao.first;
				TVar matrixval = it->second;
				z(r,ic) += (alpha*matrixval)*x.GetVal(c,ic);
			}
		}
    }
}
template<class TVar>
void TPZVerySparseMatrix<TVar>::Write(TPZStream &buf, int withclassid) const
{
	buf.Write(&this->fCol, 1);
	buf.Write(&this->fDecomposed, 1);
	buf.Write(&this->fDefPositive, 1);
	buf.Write(&this->fRow, 1);
	WriteMap(buf, withclassid, this->fExtraSparseData);
	
}
template<class TVar>
void TPZVerySparseMatrix<TVar>::WriteMap(TPZStream &buf, int withclassid, const std::map<std::pair<int64_t, int64_t>, TVar> & TheMap) const
{
	int mapsz = TheMap.size();
	buf.Write(&mapsz, 1);
	typename std::map<std::pair<int64_t, int64_t>, TVar>::const_iterator it;
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
	buf.Read(&this->fCol, 1);
	buf.Read(&this->fDecomposed, 1);
	buf.Read(&this->fDefPositive, 1);
	buf.Read(&this->fRow, 1);
	ReadMap(buf, context, this->fExtraSparseData);
	
}
template<class TVar>
void TPZVerySparseMatrix<TVar>::ReadMap(TPZStream &buf, void *context, std::map<std::pair<int64_t, int64_t>, TVar> & TheMap)
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
		std::pair<int64_t, int64_t> item(ii, jj);
		buf.Read(&value, 1);
		std::pair<std::pair<int64_t, int64_t>, TVar > fullitem(item, value);
		TheMap.insert(fullitem);
	}
}

template class TPZVerySparseMatrix<float>;
template class TPZVerySparseMatrix<double>;
template class TPZVerySparseMatrix<long double>;

template class TPZVerySparseMatrix<std::complex<float> >;
template class TPZVerySparseMatrix<std::complex<double> >;
template class TPZVerySparseMatrix<std::complex<long double> >;


