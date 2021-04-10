/**
 * @file
 * @brief Contains declaration of the abstract TPZStream class. TPZStream defines the interface for saving and reading data,\n
 * TPZFileStream implements reading from and writing to an ascii file.\n
 * TPZBFileStream implements reading from and writing to a binary file.\n
 * Finally, TPZBufferedStream implements reading from and writing to a buffer.
 */
#ifndef TPZSTREAM_H
#define TPZSTREAM_H

#include "pzreal.h"          // for REAL, TPZFlopCounter, is_arithmetic_pz
#include <stddef.h>          // for NULL
#include <complex>           // for complex
#include <string>           // for string
#include <map>               // for map
#include <set>               // for set
#include <type_traits>       // for is_same
#include <vector>            // for vector
#include "pzmanvector.h"     // for TPZManVector
#include "pzvec.h"           // for TPZVec
#include "tpzautopointer.h"  // for TPZAutoPointer
#include <inttypes.h>

#include "TPZPersistenceManager.h"

template <class T> class Fad;
template <int Num, class T> class TFad;

static uint64_t fCurrentVersion = 1; //TODO:AQUIFRANTake this away

#define TPZostream std::ostream

template<class T, int EXP>
class TPZChunkVector;

template <class T, int EXP>
class TPZAdmChunkVector;

/**
 * @ingroup save
 * @brief Defines the interface for saving and reading data. \ref save "Persistency"
 */

/**
 In fact, this class could use the facilities of the stream class of the std library
 This class is a subset of the functionality of the stream classes
 */
class TPZStream {
public:

    TPZStream() {
    }

    virtual ~TPZStream() {
    }

    virtual void Write(const bool val);

    virtual void Write(const int *p, int howMany = 1) = 0;

    virtual void Write(const unsigned int *p, int howMany = 1) = 0;
    
    virtual void Write(const int64_t *p, int howMany = 1) = 0;

    virtual void Write(const uint64_t *p, int howMany = 1) = 0;

#if defined WIN32 || defined __APPLE__
    /* On Windows7 64bits, long types have a 32 bits representation. We need to
     adjust this, calling the functions with 64 bits. */
    virtual void Write(const long *p, int howMany = 1);

    virtual void Write(const long unsigned int *p, int howMany = 1);
#endif

    virtual void Write(const float *p, int howMany = 1) = 0;

    virtual void Write(const double *p, int howMany = 1) = 0;

    virtual void Write(const long double *p, int howMany = 1);

    virtual void Write(const unsigned char *p, int howMany = 1) = 0;
    
    virtual void Write(const char *p, int howMany = 1) = 0;

    virtual void Write(const std::string *p, int howMany = 1);

    virtual void Write(const std::complex< float > *p, int howMany = 1) = 0;

    virtual void Write(const std::complex< double > *p, int howMany = 1) = 0;

    virtual void Write(const std::complex< long double > *p, int howMany = 1);


    virtual void Write(const TFad< 1, REAL > *p, int howMany = 1) = 0;
    
    virtual void Write(const TFad< 6, REAL > *p, int howMany = 1) = 0;
    
    virtual void Write(const TFad< 8, REAL > *p, int howMany = 1) = 0;
    
    virtual void Write(const TFad< 9, REAL > *p, int howMany = 1) = 0;
    
    virtual void Write(const TFad< 10, REAL > *p, int howMany = 1) = 0;

    virtual void Write(const TFad< 14, REAL > *p, int howMany = 1) = 0;
    
    virtual void Write(const Fad< float > *p, int howMany = 1) = 0;

    virtual void Write(const Fad< double > *p, int howMany = 1) = 0;

    virtual void Write(const Fad< long double > *p, int howMany = 1);

    virtual void Write(const Fad<std::complex< float > > *p, int howMany = 1);

    virtual void Write(const Fad<std::complex< double > >*p, int howMany = 1);

    virtual void Write(const Fad<std::complex< long double >> *p, int howMany = 1);


    void Write(const TPZFlopCounter *p, int howMany = 1);

    virtual void Read(bool &val);

    virtual void Read(int *p, int howMany = 1) = 0;

    virtual void Read(unsigned int *p, int howMany = 1) = 0;
    
    virtual void Read(int64_t *p, int howMany = 1) = 0;

    virtual void Read(uint64_t *p, int howMany = 1) = 0;

#if defined WIN32  || defined __APPLE__
    /* On Windows7 64bits, long types have a 32 bits representation. We need to
     adjust this, calling the functions with 64 bits. */
    
    virtual void Read(long *p, int howMany = 1);

    virtual void Read(long unsigned int *p, int howMany = 1);
#endif
    
    virtual void Read(float *p, int howMany = 1) = 0;

    virtual void Read(double *p, int howMany = 1) = 0;

    virtual void Read(long double *p, int howMany = 1);

    virtual void Read(unsigned char *p, int howMany = 1) = 0;
    
    virtual void Read(char *p, int howMany = 1) = 0;

    virtual void Read(std::string *p, int howMany = 1);

    virtual void Read(std::complex< float > *p, int howMany = 1) = 0;

    virtual void Read(std::complex< double > *p, int howMany = 1) = 0;

    virtual void Read(std::complex< long double > *p, int howMany = 1);


    virtual void Read(TFad< 1, REAL > *p, int howMany = 1) = 0;
    
    virtual void Read(TFad< 6, REAL > *p, int howMany = 1) = 0;
    
    virtual void Read(TFad< 8, REAL > *p, int howMany = 1) = 0;
    
    virtual void Read(TFad< 9, REAL > *p, int howMany = 1) = 0;

    virtual void Read(TFad< 10, REAL > *p, int howMany = 1) = 0;
    
    virtual void Read(TFad< 14, REAL > *p, int howMany = 1) = 0;

    virtual void Read(Fad< float > *p, int howMany = 1) = 0;

    virtual void Read(Fad< double > *p, int howMany = 1) = 0;

    virtual void Read(Fad< long double > *p, int howMany = 1);

    virtual void Read(Fad<std::complex< float >> *p, int howMany = 1);

    virtual void Read(Fad<std::complex< double >> *p, int howMany = 1);

    virtual void Read(Fad<std::complex< long double >> *p, int howMany = 1);



    void Read(TPZFlopCounter *p, int howMany = 1);

    //VECTORS AND ARRAYS
    template <class T>
    void Write(const TPZVec<T> &vec) {
      int64_t nel = vec.NElements();
      this->Write(&nel);
      if constexpr(
            std::is_same<std::string, T>::value ||
            is_arithmetic_pz<T>::value){
        if (nel)
          this->Write(&vec[0], vec.NElements());
      }
      else{
          for (auto c = 0; c < nel; c++) {
            vec[c].Write(*this, 0);
        }
      }
    }
    
    template <class T>
    void WritePointers(const TPZVec<TPZAutoPointer<T>> &vec) {
        uint64_t size = vec.NElements();
        this->Write(&size);
        for (uint64_t i = 0; i < size; ++i) {
            TPZPersistenceManager::WritePointer(vec[i].operator ->(), this);
        }
    }

    template <class T>
    void Write(const std::vector<T> &vec) {
        int nel = vec.size();
        this->Write(&nel);
        if constexpr(
            std::is_same<std::string, T>::value ||
            is_arithmetic_pz<T>::value){
            for (int c = 0; c < nel; c++)
                this->Write(&vec[c]);
        }
        else{
            for (int c = 0; c < nel; c++)
                vec[c].Write(*this, 0);
        }
    }

    void Write(const TPZVec<TPZFlopCounter> &vec) {
        int64_t nel = vec.NElements();
        this->Write(&nel);
        TPZVec<REAL> temp(nel);
        for (int iel = 0; iel < nel; iel++) {
            temp[iel] = vec[iel];
        }
        if (nel) this->Write(&temp[0], vec.NElements());
    }

    template <class T, class U> void Write(const std::map<T, U> &vec) {
        int64_t sz = vec.size();
        TPZManVector<T> keyVec(sz);
        TPZManVector<U> valVec(sz);
        int64_t count = 0;
        typename std::map<T, U>::const_iterator it;
        for (it = vec.begin(); it != vec.end(); it++) {
            keyVec[count] = it->first;
            valVec[count++] = it->second;
        }
        Write(keyVec);
        Write(valVec);
    }

    template <class T>
    void Write(const std::set<T> &vec) {
        if(std::is_same<std::string, T>::value ||
           is_arithmetic_pz<T>::value){
          int nel = vec.size();
          this->Write(&nel);
          typename std::set<T>::iterator it = vec.begin();
          while (it != vec.end()) {
            int val = *it;
            this->Write(&val);
            it++;
          }
        }
        else{
            DebugStop();
        }
    }

    template <class T>
    void Read(TPZVec<T> &vec) {
        Read(vec, NULL);
    }

    template <class T>
    void Read(TPZVec<T> &vec, void *context) {
        int64_t c, nc;
        this->Read(&nc, 1);
        vec.Resize(nc);
        if constexpr (std::is_same<std::string,T>::value ||
                      is_arithmetic_pz<T>::value){
            if (nc) this->Read(&vec[0], nc);
        }
        else{
          for (c = 0; c < nc; c++) {
            vec[c].Read(*this, context);
          }
        }
    }
    
    
    template <class T>
    void ReadPointers(TPZVec<TPZAutoPointer<T>> &vec) {
        uint64_t size;
        this->Read(&size,1);
        vec.resize(size);
        for (uint64_t i = 0; i < size; ++i) {
            vec[i] = TPZAutoPointerDynamicCast<T>(TPZPersistenceManager::GetAutoPointer(this));
        }
    }

    void Read(TPZVec<TPZFlopCounter> &vec) {
        int64_t nc;
        this->Read(&nc, 1);
        vec.Resize(nc);
        TPZVec<REAL> temp(nc);
        if (nc) this->Read(&temp[0], nc);
        for (int64_t ic = 0; ic < nc; ic++) {
            vec[ic] = temp[ic];
        }
    }

    template <class T>
    void Read(std::vector<T> &vec, void *context) {
        int c, nc;
        this->Read(&nc, 1);
        vec.resize(nc);
        if constexpr (std::is_same<std::string,T>::value ||
                      is_arithmetic_pz<T>::value){
          for (int i = 0; i < nc; i++) {
            this->Read(&vec[i]);
          }
        }
        else{
          for (c = 0; c < nc; c++) {
            vec[c].Read(*this, context);
          }
        }
    }

    template <int N>
    void Read(TPZManVector<REAL, N> &vec) {
        int64_t nc;
        this->Read(&nc, 1);
        vec.Resize(nc);
        if (nc) this->Read(&vec[0], nc);
    }

    template <int N>
    void Read(TPZManVector<int64_t, N> &vec) {
        int64_t nc;
        this->Read(&nc, 1);
        vec.Resize(nc);
        if (nc) this->Read(&vec[0], nc);
    }

    template <class T, class U> void Read(std::map<T, U> &vec) {
        TPZManVector<T> keyVec;
        this->Read(keyVec);
        TPZManVector<U> valVec;
        this->Read(valVec);
        int sz = keyVec.NElements();
#ifdef PZDEBUG
        if (sz != valVec.NElements()) {
            DebugStop();
        }
#endif
        int i;
        for (i = 0; i < sz; ++i) {
            vec[keyVec[i]] = valVec[i];
        }
    }

    template <class T>
    void Read(std::set<T> &vec) {
        if constexpr (std::is_same<std::string, T>::value || is_arithmetic_pz<T>::value){
        int nel;
        this->Read(&nel, 1);
        for (int i = 0; i < nel; i++) {
            T val;
            this->Read(&val);
            vec.insert(val);
        }
        }else{
            DebugStop();
        }
    }
    ////////////////

    template <class T>
    void WritePointers(const TPZVec<T *> &vec);
    
    template <class T>
    void ReadPointers(TPZVec<T *> &vec);

    template <class T>
    void WritePointers(const std::map<int, TPZAutoPointer<T>> &vec);
    
    template <class T>
    void ReadPointers(std::map<int, TPZAutoPointer<T> > &vec);
    
    template <class T>
    void WritePointers(const std::map<int, T *> &vec);
    
    template <class T>
    void ReadPointers(std::map<int, T *> &vec);
    
    template <class T>
    void WritePointers(const std::set<T *> &vec);
    
    template <class T>
    void ReadPointers(std::set<T *> &vec);

    template <class T, int EXP>
    void WritePointers(const TPZChunkVector<T *, EXP> &vec);
    
    template <class T, int EXP>
    void ReadPointers(TPZChunkVector<T *, EXP> &vec);

    template <class T, int EXP>
    void WritePointers(const TPZAdmChunkVector<T *, EXP> &vec);

    template <class T, int EXP>
    void ReadPointers(TPZAdmChunkVector<T *, EXP> &vec);

protected:

};

#include "TPZPersistenceManager.h"
template <class T> void TPZStream::WritePointers(const TPZVec<T *> &vec) {
    uint64_t nObjects = vec.NElements();
    this->Write(&nObjects);
    for (int64_t i = 0; i < nObjects; i++) {
        TPZPersistenceManager::WritePointer(vec[i], this);
    }
}

template <class T> void TPZStream::ReadPointers(TPZVec<T *> &vec) {
    uint64_t nObjects;
    this->Read(&nObjects);
    vec.Resize(nObjects);
    for (int64_t i = 0; i < nObjects; ++i) {
        vec[(const int64_t)(i)] = dynamic_cast<T *>(TPZPersistenceManager::GetInstance(this));
    }
}

template <class T>
void TPZStream::WritePointers(const std::map<int, TPZAutoPointer<T>> &map) {
    uint64_t nObjects = map.size();
    this->Write(&nObjects);
    typedef typename std::map<int, TPZAutoPointer<T>>::const_iterator map_it;
    map_it it;
    for (it = map.begin(); it != map.end(); it++) {
        this->Write(&(it->first));
        TPZPersistenceManager::WritePointer(it->second, this);
    }
}

template <class T>
void TPZStream::ReadPointers(std::map<int, TPZAutoPointer<T>> &map) {
    uint64_t nObjects;
    this->Read(&nObjects);
    int key;
    for (uint64_t i = 0; i < nObjects; ++i) {
        Read(&key);
        map[key] = TPZAutoPointerDynamicCast<T>(TPZPersistenceManager::GetAutoPointer(this));
    }
}

template <class T>
void TPZStream::WritePointers(const std::map<int, T *> &vec) {
    uint64_t nObjects = vec.size();
    this->Write(&nObjects);
    typedef typename std::map<int, T *>::const_iterator vec_it;
    vec_it it;
    for (it = vec.begin(); it != vec.end(); it++) {
        this->Write(&it->first);
        TPZPersistenceManager::WritePointer(it->second, this);
    }
}

template <class T> void TPZStream::ReadPointers(std::map<int, T *> &map) {
    uint64_t nObjects;
    this->Read(&nObjects);
    int key;
    for (uint64_t i = 0; i < nObjects; ++i) {
        Read(&key);
        map[key] = dynamic_cast<T *>(TPZPersistenceManager::GetInstance(this));
    }
}

template <class T> void TPZStream::WritePointers(const std::set<T *> &_set) {
    uint64_t nObjects = _set.size();
    this->Write(&nObjects);
    typedef typename std::set<T *>::const_iterator _set_it;
    _set_it it;
    while (it != _set.end()) {
        TPZPersistenceManager::WritePointer(it, this);
        it++;
    }
}

template <class T> void TPZStream::ReadPointers(std::set<T *> &_set) {
    uint64_t nObjects;
    this->Read(&nObjects);
    
    for(uint64_t i = 0; i < nObjects; ++i) {
        _set.insert(dynamic_cast<T *>(TPZPersistenceManager::GetInstance(this)));
    }
}

template <class T, int EXP>
void TPZStream::WritePointers(const TPZChunkVector<T *, EXP> &vec) {
    uint64_t nObjects = vec.NElements();
    this->Write(&nObjects);
    for (uint64_t i = 0; i < nObjects; i++) {
        TPZPersistenceManager::WritePointer(vec[i], this);
    }
}

template <class T, int EXP>
void TPZStream::ReadPointers(TPZChunkVector<T *, EXP> &vec) {
    uint64_t nObjects;
    this->Read(&nObjects);
    vec.Resize(nObjects);
    for (uint64_t i = 0; i < nObjects; ++i) {
        vec[i] = dynamic_cast<T *>(TPZPersistenceManager::GetInstance(this));
    }
}

template <class T, int EXP>
void TPZStream::WritePointers(const TPZAdmChunkVector<T *, EXP> &vec) {
    uint64_t nObjects = vec.NElements();
    this->Write(&nObjects);
    for (uint64_t i = 0; i < nObjects; i++) {
        TPZPersistenceManager::WritePointer(vec[i], this);
    }
    int compactScheme = as_integer(vec.fCompactScheme);
    this->Write(&compactScheme);
    Write(vec.fFree);
    Write(vec.fNFree);
}

template <class T, int EXP>
void TPZStream::ReadPointers(TPZAdmChunkVector<T *, EXP> &vec) {
    uint64_t nObjects;
    this->Read(&nObjects);
    vec.Resize(nObjects);
    for (uint64_t i = 0; i < nObjects; ++i) {
        vec[i] = dynamic_cast<T *>(TPZPersistenceManager::GetInstance(this));
    }
    int compactScheme;
    Read(&compactScheme);
    vec.fCompactScheme = (typename TPZAdmChunkVector<T *, EXP>::CompactScheme) compactScheme;
    Read(vec.fFree);
    Read(vec.fNFree);
}

#endif//TPZSTREAM_H
