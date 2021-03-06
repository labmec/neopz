/**
 * @file pzadmchunkthreadsafe.h
 * @brief Free store vector implementation (thread safe).
 */

#ifndef PZADMCHUNKTHREADSAFE_H
#define PZADMCHUNKTHREADSAFE_H

#include "pzadmchunk.h"
#include <mutex>

class TPZSavable;

template <class T,int EXP=10>
class TPZAdmChunkVectorThreadSafe : public TPZAdmChunkVector<T,EXP>
{
	public :
	TPZAdmChunkVectorThreadSafe<T,EXP> &operator=(const TPZAdmChunkVectorThreadSafe<T,EXP> &TPZAdmCh);
	
	TPZAdmChunkVectorThreadSafe(const TPZAdmChunkVectorThreadSafe<T,EXP> &AdmCh);
	
	TPZAdmChunkVectorThreadSafe(int numberofchunks = DEFAULTNUMBEROFCHUNKS);
	
	virtual ~TPZAdmChunkVectorThreadSafe();
	
	int AllocateNewElement();
	
	void SetFree(int index);
	
	int NFreeElements();
	
	int NElements() const;
	
	void CompactDataStructure(int type=2);
	
	int PrintFree(int i);
	
	void Resize(const int newsize);
	
	T &operator[](const int nelem) const;
	
	int FindObject(T *obj);
	
private:

	mutable std::mutex fAdmChunkVectorLock;
	
};

//--| IMPLEMENTATION |----------------------------------------------------------

template< class T , int EXP>
int TPZAdmChunkVectorThreadSafe<T,EXP>::NFreeElements() {
  std::scoped_lock<std::mutex> lck(fAdmChunkVectorLock);
	int resul = TPZAdmChunkVector<T,EXP>::NFreeElements();
	
	return resul;
}

template< class T , int EXP>
int TPZAdmChunkVectorThreadSafe<T,EXP>::PrintFree(int i) {
	std::scoped_lock<std::mutex> lck(fAdmChunkVectorLock);
	int resul = TPZAdmChunkVector<T,EXP>::PrintFree(i);
	
	return resul;
}

template< class T , int EXP>
int TPZAdmChunkVectorThreadSafe<T,EXP>::NElements() const{
	std::scoped_lock<std::mutex> lck(fAdmChunkVectorLock);
	int resul = TPZChunkVector<T,EXP>::NElements();
		
	return resul;
}

template< class T , int EXP>
TPZAdmChunkVectorThreadSafe<T,EXP>::TPZAdmChunkVectorThreadSafe(int numberofchunks) : TPZAdmChunkVector<T,EXP>(numberofchunks)
{
}

template< class T, int EXP >
TPZAdmChunkVectorThreadSafe<T,EXP>::~TPZAdmChunkVectorThreadSafe()
{
 
}

template< class T,int EXP >
int TPZAdmChunkVectorThreadSafe<T,EXP>::AllocateNewElement()
{
	std::scoped_lock<std::mutex> lck(fAdmChunkVectorLock);
	int resul = TPZAdmChunkVector<T,EXP>::AllocateNewElement();
	
	return resul;
}

template< class T,int EXP >
void TPZAdmChunkVectorThreadSafe<T,EXP>::SetFree(int index) {
	std::scoped_lock<std::mutex> lck(fAdmChunkVectorLock);
	TPZAdmChunkVector<T,EXP>::SetFree(index);
}

template< class T, int EXP >
void TPZAdmChunkVectorThreadSafe<T,EXP>::CompactDataStructure(int type) {
	std::scoped_lock<std::mutex> lck(fAdmChunkVectorLock);
	TPZAdmChunkVector<T,EXP>::CompactDataStructure(type);
}

template < class T,int EXP >
TPZAdmChunkVectorThreadSafe<T,EXP>::TPZAdmChunkVectorThreadSafe(const TPZAdmChunkVectorThreadSafe<T,EXP> &AdmCh) : TPZAdmChunkVector<T,EXP>(AdmCh)
{
}

template < class T , int EXP>
TPZAdmChunkVectorThreadSafe<T,EXP> & TPZAdmChunkVectorThreadSafe<T,EXP>::operator=(const TPZAdmChunkVectorThreadSafe<T,EXP> &AdmCh)
{
	std::scoped_lock<std::mutex> lck(fAdmChunkVectorLock);
	TPZAdmChunkVectorThreadSafe<T,EXP> &resul = (TPZAdmChunkVectorThreadSafe<T,EXP> &) TPZAdmChunkVector<T,EXP>::operator=(AdmCh);
	return resul;
}

template< class T, int EXP >
void TPZAdmChunkVectorThreadSafe<T,EXP>::Resize(const int newsize) {
	std::scoped_lock<std::mutex> lck(fAdmChunkVectorLock);
	TPZAdmChunkVector<T,EXP>::Resize(newsize);
}

template< class T, int EXP >
inline T &TPZAdmChunkVectorThreadSafe<T,EXP>::operator[](const int nelem) const {
	std::scoped_lock<std::mutex> lck(fAdmChunkVectorLock);
	T &resul = TPZChunkVector<T,EXP>::operator[](nelem);
	
	return resul;
}

template<class T, int EXP>
int TPZAdmChunkVectorThreadSafe<T,EXP>::FindObject(T *obj) {
	std::scoped_lock<std::mutex> lck(fAdmChunkVectorLock);
	int resul = TPZChunkVector<T,EXP>::FindObject(obj);
	return resul;
}


#endif // PZADMCHUNKTHREADSAFE_H

//--| PZ |----------------------------------------------------------------------
