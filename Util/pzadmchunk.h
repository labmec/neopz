/**
 * @file pzadmchunk.h
 * @brief Free store vector implementation.
 */

#ifndef PZADMCHUNK_H
#define PZADMCHUNK_H

#include "TPZChunkVector.h"
#include "pzstack.h"
#include "pzerror.h"

class TPZSavable;

/**
 * @ingroup util
 * @brief Implements a chunk vector with free store administration. \ref util "Utility"
 */

/** 
 * An object of this class allows the user to request a new object of
 * the type administered and allows the user to flag given elements as
 * unused.
 */
template <class T, int EXP = 10 >
class TPZAdmChunkVector : public TPZChunkVector<T, EXP> {
public:

    enum CompactScheme {
        NEVER = 0, NOW = 1, ALWAYS = 2
    };

    /** 
     * @brief Assignment operator. 
     * @param TPZAdmCh Vector which will be duplicated.
     */
    /**
     * It will copy the objects from TPZAdmCh will call the
     * assignment operator on all objects (also the freed objects).
     */
    TPZAdmChunkVector<T, EXP> & operator=(const TPZAdmChunkVector<T, EXP> &TPZAdmCh);

    /**
     * @brief Copy constructor.
     * @param AdmCh Object whose elements will be copied.
     */
    TPZAdmChunkVector(const TPZAdmChunkVector<T, EXP> &AdmCh);

    /**
     * Constructor with indication of the initial size of the chunk
     * allocation vector and the size of the chunks these sizes
     * cannot be modified during the lifecycle of the object.
     */
    /**
     * @brief Constructor.
     * @param numberofchunks Indicates how large the initial chunk
     * vector will be.
     */
    TPZAdmChunkVector(int numberofchunks = DEFAULTNUMBEROFCHUNKS);

    /** @brief Destructor */
    virtual ~TPZAdmChunkVector();

    /**
     * @brief Makes more room for new elements.
     * @return The index of a free element.
     */
    /** 
     * This method will search the list of free locations to return
     * the next free index \n in case there are no free indexes, this
     * method will increase the size of the chunk vector \n and returns
     * the allocated element.
     */
    int AllocateNewElement();

    /** 
     * @brief Indicate an element as free.
     * @note The object does not verify whether an element has been freed several times.
     * @param index The index of the element being put on the free stack.
     */
    void SetFree(int index);

    /**
     * @brief Access method to return the number of free elements.
     * @return Number of free elements.
     */
    inline int64_t NFreeElements() {
        return fFree.NElements();
    }

    /**
     * @brief Sets the method to compact the data structure based on the
     * @param type Type of compacting scheme to be used.
     */
    /**
     * parameter type:
     * \li \f$ when = 0 \f$ : never compact the data structure;
     * \li \f$ when = 1 \f$ : compact the data structure now;
     * \li \f$ when = 2 \f$ : compact the data structure always (default).
     */
    void CompactDataStructure(CompactScheme type = CompactScheme::ALWAYS );

    /** @brief Print index i into the fFree vector. */
    inline int PrintFree(int i) {
        return fFree[i];
    }

    /**
     * @brief  Increase the size of the chunk vector.
     * @param newsize Requested new size of the vector.
     */
    void Resize(const int newsize);

    int ClassId() const override {
        return Hash("TPZAdmChunkVector") ^ TPZChunkVector<T, EXP>::ClassId() << 1;
    }
    
    void Read(TPZStream& buf, void* context) override{
        uint64_t c, nc;
        buf.Read(&nc, 1);
        this->Resize(nc);
        for (c = 0; c < nc; c++){
            ReadInternal(this->operator [](c), buf, context);
        }
        int compactScheme;
        buf.Read(&compactScheme, 1);
        this->fCompactScheme = TPZAdmChunkVector<T, EXP>::CompactScheme(compactScheme);
        buf.Read(this->fFree);
        buf.Read(this->fNFree);
    }
    
    void Write(TPZStream& buf, int withclassid) const override{
        uint64_t c, nc = this->NElements();
        buf.Write(&nc);
        for (c = 0; c < nc; c++){
            WriteInternal(this->operator [](c), buf, withclassid);
        }
        int val = as_integer( this->fCompactScheme );
        buf.Write(&val);
        buf.Write(this->fFree);
        buf.Write(this->fNFree);
    }

private:

    friend class TPZStream;
    friend class TPZGeoMesh;

    /**
     * @brief Internal variable indicating the type of compacting scheme.
     * @see CompactDataStructure.
     */
    CompactScheme fCompactScheme;

    /** @brief Number of free elements within each chunk. */
    TPZManVector<int> fNFree;

    /** @brief List of indexes of freed elements. */
    TPZStack<int> fFree;
};

//--| IMPLEMENTATION |----------------------------------------------------------

template< class T, int EXP>
TPZAdmChunkVector<T, EXP>::TPZAdmChunkVector(int numberofchunks)
: TPZChunkVector<T, EXP>(numberofchunks),
fCompactScheme(NEVER), // never compact the data structure
fNFree(numberofchunks),
fFree() {
    for (int i = 0; i < numberofchunks; i++) {
        fNFree[i] = 0;
    }

    fNFree.Resize(0);
}

template< class T, int EXP >
TPZAdmChunkVector<T, EXP>::~TPZAdmChunkVector() {
}

// Return the index of a free element

template< class T, int EXP >
int TPZAdmChunkVector<T, EXP>::AllocateNewElement() {
    if (fFree.NElements() > 0) {
        int index = fFree.Pop();
        int chunk = index >> EXP;
        fNFree[chunk]--;
        return index;
    }

    Resize(this->NElements() + 1);

    return this->NElements() - 1;
}

// Indicate an element as free
template< class T, int EXP >
void TPZAdmChunkVector<T, EXP>::SetFree(int index) {
#ifndef NODEBUG
    if (index < 0) {
        PZError << "TPZAdmChunkVector::SetFree. Bad parameter index." << std::endl;
        PZError.flush();
        return;
    }
#endif

    int chunk = index >> EXP;

    fNFree[chunk]++;
    fFree.Push(index);

    if (fCompactScheme == ALWAYS) {
        CompactDataStructure(ALWAYS);
}
}

// Let to compact the data structure.
// If type=0 never compact, (type=1 let to compact now),(type=2  let to compact always)

template< class T, int EXP >
void TPZAdmChunkVector<T, EXP>::CompactDataStructure(CompactScheme type) {
    switch (type) {
        case NEVER:
            fCompactScheme = NEVER;
            return;
        case ALWAYS:
            fCompactScheme = ALWAYS;
        case NOW:
        {
            int chunksize = 1 << EXP;
            int nchunksused = 0;
            if (this->NElements()) nchunksused = ((this->NElements() - 1) >> EXP) + 1;
            int i = nchunksused - 1;
            int maxfree = this->NElements()-((nchunksused - 1) << EXP);

            if (i >= 0 && this->fVec[i] && fNFree[i] == maxfree) {
                Resize(chunksize * i);
                --i;
                while (i >= 0 && this->fVec[i] && fNFree[i] == chunksize) {
                    Resize(chunksize * i);
                    --i;
                }
            }
            this->fVec.Shrink();
            fNFree.Shrink();
            fFree.Shrink();
            break;
        }
        default:
            PZError << "TPZAdmChunkVector::CompactDataStructure. Bad parameter type."
                    << std::endl;
            PZError.flush();
    }
}

template < class T, int EXP >
TPZAdmChunkVector<T, EXP>::TPZAdmChunkVector(const TPZAdmChunkVector<T, EXP> &AdmCh) :
TPZChunkVector<T, EXP>(AdmCh), fCompactScheme(AdmCh.fCompactScheme),
fNFree(AdmCh.fNFree), fFree(AdmCh.fFree) {
    // NOTHING TO DO HERE!
}

template < class T, int EXP>
TPZAdmChunkVector<T, EXP> & TPZAdmChunkVector<T, EXP>::operator=(
        const TPZAdmChunkVector<T, EXP> &AdmCh) {
    if (this == &AdmCh)
        return *this;

    TPZChunkVector<T, EXP>::operator=(AdmCh);

    fFree = AdmCh.fFree;
    fNFree = AdmCh.fNFree;
    fCompactScheme = AdmCh.fCompactScheme;

    return *this;
}

template< class T, int EXP >
void TPZAdmChunkVector<T, EXP>::Resize(const int newsize) {
#ifndef NODEBUG
    if (newsize < 0) {
        PZError << "TPZAdmChunkVector::Resize. Bad parameter newsize." << std::endl;
        PZError.flush();
        return;
    }
#endif

    TPZChunkVector<T, EXP>::Resize(newsize);

    //   int sizechunk = 1 << EXP;
    int nchunks = fNFree.NElements();
    int chunksneeded = this->fVec.NElements(); // equivalent to newsize>>fExponent??

    fNFree.Resize(chunksneeded);

    for (int i = nchunks; i < chunksneeded; i++) {
        fNFree[i] = 0;
    }

    if (chunksneeded > nchunks) return;

    // delete all free indexes which are above the new size
    // update the number of free elements of the last chunk
    TPZStack<int> temp(fFree);
    temp.Resize(0);

    while (fFree.NElements() > 0) {
        int index = fFree.Pop();
        if (index < newsize)
            temp.Push(index);
        else {
            int chunk = index >> EXP;

            if (chunk == chunksneeded - 1)
                fNFree[chunksneeded - 1]--;
        }
    }

    fFree = temp;
}

#endif // PZADMCHUNK_H

//--| PZ |----------------------------------------------------------------------
