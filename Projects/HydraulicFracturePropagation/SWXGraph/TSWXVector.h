
#ifndef TSWXVectorH
#define TSWXVectorH

#include <vector>
#include <iostream>
#include <math.h>
#include <map>

#include "pzerror.h"
#include "pzreal.h"

namespace swx{

template< class T >
class vector{

private:

    std::vector<T> fStore;

public:

  vector():fStore(){
    ///nothing
  }

  vector(unsigned int size):fStore(size){
    ///nothing
  }

  vector(const int size, const T& copy):fStore(size,copy){
#ifdef PZDEBUG
    if( size < 0 ) DebugStop();
#endif
  }

  vector(const vector<T> &copy):fStore(){
    this->operator=(copy);
  }

  virtual ~vector(){
    ///nothing
  }

  vector<T> &operator=(const vector<T> &copy){
    const int n = copy.size();
	this->resize( n );
    for(unsigned int i = 0; i < n; i++){
      this->operator[](i) = copy[i];
    }
    return *this;
  }

  bool operator==(const vector<T> &other) const;

  std::vector<T> & GetStdVector(){
	return fStore;
  }

  const std::vector<T> & GetStdVector() const{
    return fStore;
  }

  void Fill(const T& a){
    for(unsigned int i = 0; i < fStore.size(); i++) fStore[i] = a;
  }

  void clear(){
    fStore.clear();
  }

  const T & operator[]( int index ) const{
#ifdef PZDEBUG
    unsigned int size = fStore.size();
    if( index < 0 || index >= size ){
	  DebugStop();
    }
#endif
    return fStore[ index ];
  }///method

  T & operator[]( int index ){
#ifdef PZDEBUG
		unsigned int size = fStore.size();
		if( index < 0 || index >= size ){
			DebugStop();
		}
#endif
    return fStore[index];
  }///method

  unsigned int size() const{
    return fStore.size();
  }

	void resize(int newsize){
#ifdef PZDEBUG
	if( newsize < 0 ) DebugStop();
#endif
	fStore.resize(newsize);
	}

	void erase(int index)
	{
		int size = this->size();
#ifdef PZDEBUG
		if( index < 0 || index >= size ){
			DebugStop();
		}
#endif
		for(unsigned p = index+1; p < size; p++)
		{
			fStore[p-1] = fStore[p];
		}
		fStore.resize(size-1);
	}

  void push_back(const T &a){
  	fStore.push_back(a);
  }

  void Write(std::ostream &file) const;
  void Read(std::istream &file);
#if defined(SWX_BUILDER_2010) || defined (SWX_BUILDER_XE2)
  void Write(_di_IXMLNode &myNode) const;
  void Read(_di_IXMLNode &myNode);
#endif
  void insert(int pos, const T &val){
    typename std::vector<T>::iterator w = fStore.begin();
    w+=pos;
    fStore.insert(w,val);
  }///void

};///class

template< >
inline void vector<double>::Write(std::ostream &file) const{
  file << this->size() << "\t";
  for(unsigned int i = 0; i < this->size(); i++){
    file << this->operator[](i) << "\t";
  }
  file << "\n";
}///void

template< >
inline void vector<int>::Write(std::ostream &file) const{
  file << this->size() << "\t";
  for(unsigned int i = 0; i < this->size(); i++){
    file << this->operator[](i) << "\t";
  }
  file << "\n";
}///void
#if defined(SWX_BUILDER_2010) || defined(SWX_BUILDER_XE2)
template< >
inline void vector< System::UnicodeString >::Write(std::ostream &file) const{
  file << this->size() << "\t";
  for(unsigned int i = 0; i < this->size(); i++){
		std::string s;
#ifdef SWX_BUILDER_2010
		s = System::UnicodeString(this->operator[](i)).t_str();
#endif
#ifdef SWX_BUILDER_XE2
    s = swx::wstring2string(this->operator[](i).c_str());
#endif
    file << s << "\t";
  }
  file << "\n";
}///void
#endif
template< class T >
inline void vector< T >::Write(std::ostream &file) const{
  DebugStop();
}///void

#if defined(SWX_BUILDER_2010) || defined (SWX_BUILDER_XE2)
template< class T >
inline void vector< T >::Write(_di_IXMLNode &myNode) const{
  std::stringstream writeSSTR;
  this->Write(writeSSTR);
#ifdef SWX_BUILDER_2010
  myNode->Text = writeSSTR.str().c_str();
#endif
#ifdef SWX_BUILDER_XE2
  myNode->Text = swx::string2wstring(writeSSTR.str()).c_str();
#endif
}
#endif
template< >
inline void vector<double>::Read(std::istream &file){
  unsigned int n;
  file >> n;
  this->resize(n);
  for(unsigned int i = 0; i < n; i++){
    file >> this->operator[](i);
  }
}///void

template< >
inline void vector<int>::Read(std::istream &file){
  unsigned int n;
  file >> n;
  this->resize(n);
  for(unsigned int i = 0; i < n; i++){
    file >> this->operator[](i);
  }
}///void
#if defined(SWX_BUILDER_2010) || defined(SWX_BUILDER_XE2)
template< >
inline void vector< System::UnicodeString >::Read(std::istream &file){
  unsigned int n;
  file >> n;
  this->resize(n);
  for(unsigned int i = 0; i < n; i++){
    std::string s;
    file >> s;
#ifdef SWX_BUILDER_2010
    this->operator[](i) = s.c_str();
#endif
#ifdef SWX_BUILDER_XE2
    this->operator[](i) = swx::string2wstring(s).c_str();
#endif
  }
}///void
#endif
template< class T >
inline void vector< T >::Read(std::istream &file){
  DebugStop();
}///void

#if defined(SWX_BUILDER_2010) || defined (SWX_BUILDER_XE2)
template< class T >
inline void vector< T >::Read(_di_IXMLNode &myNode){
  std::stringstream readSSTR;
#ifdef SWX_BUILDER_2010
  readSSTR << myNode->Text.t_str();
#endif
#ifdef SWX_BUILDER_XE2
  readSSTR << swx::wstring2string(myNode->Text.c_str());
#endif
  this->Read(readSSTR);
}
#endif

template< >
inline bool vector<double>::operator==(const vector<double> &other) const{
  if(this->size() != other.size()) return false;
  for(unsigned int i = 0; i < this->size(); i++){
	if( fabs(this->operator[](i) - other[i]) > 1e-16 ) return false;
  }
  return true;
}

template< class T >
inline bool vector< T >::operator==(const vector< T > &other) const{
  if(this->size() != other.size()) return false;
  for(unsigned int i = 0; i < this->size(); i++){
	if( !(this->operator[](i) == other[i]) ) return false;
  }
  return true;
}

template< class T >
inline void WriteVec(std::ostream &file, const swx::vector< T > &vec ){
  file << vec.size() << "\n";
  for(unsigned int i = 0; i < vec.size(); i++){
    vec[i].Write(file);
  }
}

template< class T >
inline void ReadVec(std::istream &file, swx::vector< T > &vec ){
  int n;
  file >> n;
  vec.resize( n );
  for(unsigned int i = 0; i < vec.size(); i++){
    vec[i].Read(file);
  }
}

};///namespace

//---------------------------------------------------------------------------
#endif


