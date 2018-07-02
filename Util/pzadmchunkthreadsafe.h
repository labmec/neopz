/**
 * @file pzadmchunkthreadsafe.h
 * @brief Free store vector implementation (thread safe).
 */

#ifndef PZADMCHUNKTHREADSAFE_H
#define PZADMCHUNKTHREADSAFE_H

#include "pzadmchunk.h"
#include <pthread.h>

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

	mutable pthread_mutex_t fAdmChunkVectorLock;
	
};

//--| IMPLEMENTATION |----------------------------------------------------------

template< class T , int EXP>
int TPZAdmChunkVectorThreadSafe<T,EXP>::NFreeElements() {
	pthread_mutex_lock(&fAdmChunkVectorLock);
	int resul = TPZAdmChunkVector<T,EXP>::NFreeElements();
	pthread_mutex_unlock(&fAdmChunkVectorLock);
	
	return resul;
}

template< class T , int EXP>
int TPZAdmChunkVectorThreadSafe<T,EXP>::PrintFree(int i) {
	pthread_mutex_lock(&fAdmChunkVectorLock);
	int resul = TPZAdmChunkVector<T,EXP>::PrintFree(i);
	pthread_mutex_unlock(&fAdmChunkVectorLock);
	
	return resul;
}

template< class T , int EXP>
int TPZAdmChunkVectorThreadSafe<T,EXP>::NElements() const{
	pthread_mutex_lock(&fAdmChunkVectorLock);
	int resul = TPZChunkVector<T,EXP>::NElements();
	pthread_mutex_unlock(&fAdmChunkVectorLock);
	
	return resul;
}

template< class T , int EXP>
TPZAdmChunkVectorThreadSafe<T,EXP>::TPZAdmChunkVectorThreadSafe(int numberofchunks) : TPZAdmChunkVector<T,EXP>(numberofchunks)
{
	pthread_mutex_init(&fAdmChunkVectorLock, NULL);
}

template< class T, int EXP >
TPZAdmChunkVectorThreadSafe<T,EXP>::~TPZAdmChunkVectorThreadSafe()
{
	pthread_mutex_destroy(&fAdmChunkVectorLock);
}

template< class T,int EXP >
int TPZAdmChunkVectorThreadSafe<T,EXP>::AllocateNewElement()
{
	pthread_mutex_lock(&fAdmChunkVectorLock);
	int resul = TPZAdmChunkVector<T,EXP>::AllocateNewElement();
	pthread_mutex_unlock(&fAdmChunkVectorLock);
	
	return resul;
}

template< class T,int EXP >
void TPZAdmChunkVectorThreadSafe<T,EXP>::SetFree(int index) {
	pthread_mutex_lock(&fAdmChunkVectorLock);
	TPZAdmChunkVector<T,EXP>::SetFree(index);
	pthread_mutex_unlock(&fAdmChunkVectorLock);
}

template< class T, int EXP >
void TPZAdmChunkVectorThreadSafe<T,EXP>::CompactDataStructure(int type) {
	pthread_mutex_lock(&fAdmChunkVectorLock);
	TPZAdmChunkVector<T,EXP>::CompactDataStructure(type);
	pthread_mutex_unlock(&fAdmChunkVectorLock);
}

template < class T,int EXP >
TPZAdmChunkVectorThreadSafe<T,EXP>::TPZAdmChunkVectorThreadSafe(const TPZAdmChunkVectorThreadSafe<T,EXP> &AdmCh) : TPZAdmChunkVector<T,EXP>(AdmCh)
{
	pthread_mutex_init(&fAdmChunkVectorLock, NULL);
}

template < class T , int EXP>
TPZAdmChunkVectorThreadSafe<T,EXP> & TPZAdmChunkVectorThreadSafe<T,EXP>::operator=(const TPZAdmChunkVectorThreadSafe<T,EXP> &AdmCh)
{
	pthread_mutex_lock(&fAdmChunkVectorLock);
	TPZAdmChunkVectorThreadSafe<T,EXP> &resul = (TPZAdmChunkVectorThreadSafe<T,EXP> &) TPZAdmChunkVector<T,EXP>::operator=(AdmCh);
	pthread_mutex_unlock(&fAdmChunkVectorLock);
	return resul;
}

template< class T, int EXP >
void TPZAdmChunkVectorThreadSafe<T,EXP>::Resize(const int newsize) {
	pthread_mutex_lock(&fAdmChunkVectorLock);
	TPZAdmChunkVector<T,EXP>::Resize(newsize);
	pthread_mutex_unlock(&fAdmChunkVectorLock);
}

template< class T, int EXP >
inline T &TPZAdmChunkVectorThreadSafe<T,EXP>::operator[](const int nelem) const {
	pthread_mutex_lock(&fAdmChunkVectorLock);
	T &resul = TPZChunkVector<T,EXP>::operator[](nelem);
	pthread_mutex_unlock(&fAdmChunkVectorLock);
	
	return resul;
}

template<class T, int EXP>
int TPZAdmChunkVectorThreadSafe<T,EXP>::FindObject(T *obj) {
	pthread_mutex_lock(&fAdmChunkVectorLock);
	int resul = TPZChunkVector<T,EXP>::FindObject(obj);
	pthread_mutex_unlock(&fAdmChunkVectorLock);
	
	return resul;
}


#endif // PZADMCHUNKTHREADSAFE_H

//--| PZ |----------------------------------------------------------------------
