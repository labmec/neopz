/**
 * @file
 * @brief Contains the TPZMatWithMem class which implements the memory features.
 */

#ifndef PZMATWITHMEM_H
#define PZMATWITHMEM_H

#include "TPZMaterial.h"
#include "pzadmchunk.h"
#include <memory>

/**
 * @ingroup material
 * @brief Implements an abstract class implementing the memory features.
 */

/**
 * Materials that aim to use memory shoud be derived from this class.
 */

template <class TMEM, class TFather = TPZMaterial>
class TPZMatWithMem : public TFather {
public:

    /** @brief Default constructor */
    TPZMatWithMem();

    /** @brief Creates a material object and inserts it in the vector of material pointers of the mesh */
    /** Upon return vectorindex contains the index of the material object within the vector */
    TPZMatWithMem(int id);

    /** 
     * @brief Creates a material object based on the referred object and
     * inserts it in the vector of material pointers of the mesh
     */
    /** Upon return vectorindex contains the index of the material object within the vector */
    TPZMatWithMem(const TPZMatWithMem<TMEM, TFather> &mat);

    virtual ~TPZMatWithMem();

    /** @brief Returns the name of the material */
    virtual std::string Name() {
        return "TPZMatWithMem< >";
    }

    /** @brief Prints out the data associated with the material */
    virtual void PrintMem(std::ostream &out = std::cout, const int memory = 0);

    /** @brief Prints out the data associated with the material */
    virtual void Print(std::ostream &out);

    virtual TMEM &MemItem(const int i) const;

public:
    
    /** @brief Unique identifier for serialization purposes */
    virtual int ClassId() const;

    virtual void Write(TPZStream &buf, int withclassid) const;

    virtual void Read(TPZStream &buf, void *context);

    std::shared_ptr<TPZAdmChunkVector<TMEM>> & GetMemory();

    void SetMemory(std::shared_ptr<TPZAdmChunkVector<TMEM>> & memory);

    /**
     * @brief Pushes a new entry in the context of materials with memory
     * @return Returning its index at the internal storage stack
     */
    /** To be implemented only in the proper materials. */
    virtual int PushMemItem(int sourceIndex = -1);

    /** @brief Frees an entry in the material with memory internal history storage */
    virtual void FreeMemItem(int index);

    /** @ Reset the memory index to its default value */
    void ResetMemItem(int index) {
        this->MemItem(index) = fDefaultMem;
    }

    /// Reset all memory items
    void ResetMemory() {
        int nmem = fMemory->NElements();
        for (unsigned int im = 0; im < nmem; im++) {
            ResetMemItem(im);
        }
    }

    /** @brief Sets the default memory settings for initialization */
    virtual void SetDefaultMem(TMEM & defaultMem);

    /** @brief Sets/Unsets the internal memory data to be updated in the next assemble/contribute call */
    virtual void SetUpdateMem(bool update = 1);

protected:

    /** @brief Shared pointer to material memory */
    std::shared_ptr<TPZAdmChunkVector<TMEM>> fMemory;

    /** @brief Default memory settings */
    TMEM fDefaultMem;

    /** @brief Flag to indicate whether the memory data are to be updated in an assemble loop */
    bool fUpdateMem;
};

template <class TMEM, class TFather>
TPZMatWithMem<TMEM, TFather>::TPZMatWithMem() :
TPZRegisterClassId(&TPZMatWithMem::ClassId),
TFather(),
fMemory(new TPZAdmChunkVector<TMEM>()), fDefaultMem(), fUpdateMem(0) {
}

template <class TMEM, class TFather>
TPZMatWithMem<TMEM, TFather>::TPZMatWithMem(int id) :
TPZRegisterClassId(&TPZMatWithMem::ClassId),
TFather(id),
fMemory(new TPZAdmChunkVector<TMEM>()), fDefaultMem(), fUpdateMem(0) {
}

template <class TMEM, class TFather>
TPZMatWithMem<TMEM, TFather>::TPZMatWithMem(const TPZMatWithMem<TMEM, TFather> &mat) :
TPZRegisterClassId(&TPZMatWithMem::ClassId),
TFather(mat),
fMemory(new TPZAdmChunkVector<TMEM>()), fDefaultMem(mat.fDefaultMem), fUpdateMem(mat.fUpdateMem) {
    *fMemory = *mat.fMemory;
}

template <class TMEM, class TFather>
TPZMatWithMem<TMEM, TFather>::~TPZMatWithMem() {
}

template <class TMEM, class TFather>
void TPZMatWithMem<TMEM, TFather>::Print(std::ostream &out) {
    out << __PRETTY_FUNCTION__ << std::endl;
    TFather::Print(out);

    out << "\nfDefaultMem = \n" << fDefaultMem;
    out << "\nfUpdateMem = " << fUpdateMem;
    int size = fMemory->NElements();
    out << "\nfMemory with " << size << " elements";
    for (int i = 0; i < size; i++) {
        out << "\nfMemory element : " << i << std::endl;
        this->MemItem(i).Print(out);
    }

}

template <class TMEM, class TFather>
void TPZMatWithMem<TMEM, TFather>::PrintMem(std::ostream &out, const int memory) {
    int size = fMemory->NElements();
    if (memory >= 0 && memory < size) {
        out << "fMemory element : " << memory << std::endl;
        this->MemItem(memory).Print(out);
    } else {
        out << "Memory index out of range : memory " << memory << " no elements " << size << std::endl;
    }
}

template <class TMEM, class TFather>
TMEM &TPZMatWithMem<TMEM, TFather>::MemItem(const int i) const {
    return fMemory.get()->operator [](i);
}

template <class TMEM, class TFather>
int TPZMatWithMem<TMEM, TFather>::ClassId() const {
    return Hash("TPZMatWithMem") ^ TFather::ClassId() << 1 ^ TMEM().ClassId() << 2;
}

template <class TMEM, class TFather>
void TPZMatWithMem<TMEM, TFather>::Write(TPZStream &buf, int withclassid) const {
    TFather::Write(buf, withclassid);
    int updatemem = fUpdateMem;
    buf.Write(&updatemem);
    fDefaultMem.Write(buf, 0);
    TPZPersistenceManager::WritePointer(fMemory.get(), &buf);
}

template <class TMEM, class TFather>
void TPZMatWithMem<TMEM, TFather>::Read(TPZStream &buf, void *context) {
    TFather::Read(buf, context);
    int updatemem;
    buf.Read(&updatemem);
    if (updatemem) {
        fUpdateMem = true;
    } else {
        fUpdateMem = false;
    }
    fDefaultMem.Read(buf, 0);
    fMemory = std::dynamic_pointer_cast<TPZAdmChunkVector<TMEM> >(TPZPersistenceManager::GetSharedPointer(&buf));

}

template <class TMEM, class TFather>
std::shared_ptr<TPZAdmChunkVector<TMEM>> & TPZMatWithMem<TMEM, TFather>::GetMemory() {
    return fMemory;
}

template <class TMEM, class TFather>
void TPZMatWithMem<TMEM, TFather>::SetMemory(std::shared_ptr<TPZAdmChunkVector<TMEM>> & memory) {
    fMemory = memory;
}

template <class TMEM, class TFather>
int TPZMatWithMem<TMEM, TFather>::PushMemItem(int sourceIndex) {
    int index = fMemory->AllocateNewElement();
    if (sourceIndex < 0) {
        this->ResetMemItem(index);
    } else {
        this->MemItem(index) = this->MemItem(sourceIndex);
    }
    return index;
}

template <class TMEM, class TFather>
void TPZMatWithMem<TMEM, TFather>::FreeMemItem(int index) {
    fMemory->SetFree(index);
}

template <class TMEM, class TFather>
void TPZMatWithMem<TMEM, TFather>::SetDefaultMem(TMEM & defaultMem) {
    fDefaultMem = defaultMem;
}

template <class TMEM, class TFather>
void TPZMatWithMem<TMEM, TFather>::SetUpdateMem(bool update) {
    fUpdateMem = update;
}

#endif
