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
int TPZVerySparseMatrix<TVar>::PutVal(const long row,const long col, const TVar &val)
{
    if (row < 0 || col < 0 || row >this->fRow || col >this->fCol)
    {
        cout<< "ERRO! em TPZVerySparseMatrix::PutVal: The row i or column j are incompatible with the rows or columns of a matrix"<< endl;
        return -1;
	}
	
    pair <long,long> position(row,col);
	typename std::map <std::pair<long, long>, TVar>::iterator it = this->fExtraSparseData.find(position);
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
TPZVerySparseMatrix<TVar>::TPZVerySparseMatrix(const TPZFMatrix<TVar> &cp) : TPZMatrix<TVar>(cp)
{
	for(long i=0; i<this->fRow; i++)
	{
		for(long j=0; j<this->fCol; j++)
		{
			TVar a = cp.GetVal(i,j);
			if(!IsZero(a)) PutVal(i,j,a);
		}
	}
}

template<class TVar>
void TPZVerySparseMatrix<TVar>::Simetrize() {
  
  long rows = this->Rows();
  long cols = this->Cols();
  
  // TODO: introduce error handling mechanism  
  if ( rows != cols) {
      cerr << "Error: Simetrize only work for square matrices";
      return;
  }
  
  typename std::map <std::pair<long, long>, TVar>::iterator it = this->fExtraSparseData.begin();
  typename std::map <std::pair<long, long>, TVar>::iterator end = this->fExtraSparseData.end();
  typename std::map <std::pair<long, long>, TVar>::iterator next;
  
  std::list< std::pair< std::pair<long, long>, TVar > > temp;
  
  for(; it != end; it=next) {
    const std::pair<long, long>& key = it->first;
    next = it;
    next++;
    
    if(key.first < key.second) {
      temp.push_back( std::pair< std::pair<long, long>, TVar > (std::pair<long, long>(key.second, key.first), it->second));
    }
    else if (key.first > key.second) {
      this->fExtraSparseData.erase(it);
    }
  }
  
  typename std::list< std::pair< std::pair<long, long>, TVar > >::iterator at = temp.begin();
  typename std::list< std::pair< std::pair<long, long>, TVar > >::iterator atEnd = temp.end();
  
  	for(; at != atEnd; at++) {
		this->fExtraSparseData [ at->first ] = at->second;
	}

}

template<class TVar>
void TPZVerySparseMatrix<TVar>::Transpose(TPZVerySparseMatrix<TVar> *T) const {
  long rows = this->Rows();
  long cols = this->Cols();
  T->Resize( cols, rows );
  T->fExtraSparseData.clear();
  typename std::map <std::pair<long, long>, TVar>::const_iterator it = this->fExtraSparseData.begin();
  typename std::map <std::pair<long, long>, TVar>::const_iterator end = this->fExtraSparseData.end();

  for (; it != end; it++) {
    const std::pair<long, long>& key = it->first;
    T->fExtraSparseData[std::pair<long,long>(key.second, key.first)] = it->second;
  }
}

template<class TVar>
const TVar & TPZVerySparseMatrix<TVar>::GetVal(const long row, const long col) const
{
    if (row < 0 || col < 0 || row >this->fRow || col >this->fCol)
    {
        cout<< "ERRO! em TPZVerySparseMatrix::GetVal: The row i or column j are incompatible with the rows or columns of a matrix"<< endl;
		return this->gZero;
    }
	
    pair <long,long> position(row,col);
    typename map<pair<long,long>, TVar>::const_iterator it;
    it = fExtraSparseData.find(position);
	
    if (it == fExtraSparseData.end() )
    {
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
    
    long xcols = x.Cols();
    long ic, c, r;
    this->PrepareZ(y,z,beta,opt,stride);
    TVar val = 0.;
	
    for (ic = 0; ic < xcols; ic++)
    {
        if(!opt) 
        {
			typename map< pair<long,long>, TVar>::const_iterator it;
			
			for(it = fExtraSparseData.begin(); it!= fExtraSparseData.end(); it++)
			{
				pair <long, long> position(it->first);
				c = position.second;
				r = position.first;
				TVar matrixval = it->second;
				
				val = z(r*stride,ic) + alpha * matrixval * x.GetVal(c*stride,ic);
				z.PutVal(r*stride,ic,val);
			}
        } else 
		{
			typename map<pair<long,long>, TVar>::const_iterator it;
			
			for(it = fExtraSparseData.begin(); it != fExtraSparseData.end(); it++)
			{
				pair <long, long> posicao(it->first);
				r = posicao.second;
				c = posicao.first;
				TVar matrixval = it->second;
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
void TPZVerySparseMatrix<TVar>::WriteMap(TPZStream &buf, int withclassid, std::map<std::pair<long, long>, TVar> & TheMap)
{
	int mapsz = TheMap.size();
	buf.Write(&mapsz, 1);
	typename std::map<std::pair<long, long>, TVar>::iterator it;
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
void TPZVerySparseMatrix<TVar>::ReadMap(TPZStream &buf, void *context, std::map<std::pair<long, long>, TVar> & TheMap)
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
		std::pair<long, long> item(ii, jj);
		buf.Read(&value, 1);
		std::pair<std::pair<long, long>, TVar > fullitem(item, value);
		TheMap.insert(fullitem);
	}
}

template class TPZVerySparseMatrix<float>;
template class TPZVerySparseMatrix<double>;
template class TPZVerySparseMatrix<long double>;

template class TPZVerySparseMatrix<std::complex<float> >;
template class TPZVerySparseMatrix<std::complex<double> >;
template class TPZVerySparseMatrix<std::complex<long double> >;


