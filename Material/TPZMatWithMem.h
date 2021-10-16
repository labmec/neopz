#ifndef TPZMATWITHMEM_H
#define TPZMATWITHMEM_H

#include "pzadmchunk.h"
#include <memory>

//! Type-agnostic base class for TPZMatWithMem
class TPZMatWithMemBase : public virtual TPZSavable{
public:
    //! Default constructor
    TPZMatWithMemBase() = default;

    /** @brief Prints out the data associated with the material */
    virtual void PrintMem(std::ostream &out = std::cout, const int memory = 0) const= 0;

    /** @brief Prints out the data associated with the material */
    virtual void Print(std::ostream &out) const= 0;
    /**
     * @brief Pushes a new entry in the context of materials with memory
     * @return Returning its index at the internal storage stack
     */

    virtual int PushMemItem(int sourceIndex = -1) = 0;

    /** @brief Frees an entry in the material with memory internal history storage */
    virtual void FreeMemItem(int index) = 0;

    /** @ Reset the memory index to its default value */
    virtual void ResetMemItem(int index) = 0;       

    /// Reset all memory items
    virtual void ResetMemory() = 0;

    /** @brief Sets/Unsets the internal memory data to be updated in the next assemble/contribute call */
    virtual void SetUpdateMem(bool update = 1) = 0;
    
    /** @brief Gets the internal memory data to be updated in the next assemble/contribute call */
    virtual bool GetUpdateMem() = 0 ;
};


template<class TMem>
class TPZMatWithMemBC;
/**
   @ingroup singlespaceinterface combinedspacesinterface
   @brief Defines an interface for materials with memory in the integration points.
 */
template <class TMem>
class TPZMatWithMem : public TPZMatWithMemBase {
public:
    // this is type alias
    // https://en.cppreference.com/w/cpp/language/type_alias
    // from now on we can use TPZMatCombinedSpacesT<TVar>::TInterfaceBC as a type
    // this will be used in CreateBC
    using TInterfaceBC = TPZMatWithMemBC<TMem>;

    //! Default constructor
    TPZMatWithMem();
    //! Destructor
    ~TPZMatWithMem() = default;
    //! Copy constructor
    TPZMatWithMem(const TPZMatWithMem &mat);
    //! Move constructor
    TPZMatWithMem(TPZMatWithMem &&mat) = default;
    //! Copy-assignment operator
    TPZMatWithMem& operator=(const TPZMatWithMem &mat);
    //! Move-assignment operator
    TPZMatWithMem& operator=(TPZMatWithMem &&mat) = default;
    //!@name Memory
    /** @{*/
    /** @brief Prints out the data associated with the material */
    void PrintMem(std::ostream &out = std::cout, const int memory = 0) const override;
    //@}
    /** @brief Prints out the data associated with the material */
    void Print(std::ostream &out) const override;
    //!@name Memory
    /** @{*/
    virtual TMem &MemItem(const int i) const;
    /** @}*/
    /** @brief Unique identifier for serialization purposes */
    int ClassId() const override;

    void Write(TPZStream &buf, int withclassid) const override;

    void Read(TPZStream &buf, void *context) override;
    //!@name Memory
    /** @{*/
    std::shared_ptr<TPZAdmChunkVector<TMem>> & GetMemory();

    void SetMemory(std::shared_ptr<TPZAdmChunkVector<TMem>> & memory);

    /**
     * @brief Pushes a new entry in the context of materials with memory
     * @return Returning its index at the internal storage stack
     */
    int PushMemItem(int sourceIndex = -1) override;

    /** @brief Frees an entry in the material with memory internal history storage */
    void FreeMemItem(int index) override;

    /** @brief Reset the memory index to its default value */
    void ResetMemItem(int index) override{
        this->MemItem(index) = fDefaultMem;
    }

    /// Reset all memory items
    void ResetMemory() override{
        int nmem = fMemory->NElements();
        for (unsigned int im = 0; im < nmem; im++) {
            ResetMemItem(im);
        }
    }

    /** @brief Sets the default memory settings for initialization */
    virtual void SetDefaultMem(TMem & defaultMem);

    /// Access the default memory settings for initislization
    virtual TMem &GetDefaultMemory()
    {
        return fDefaultMem;
    }

    /** @brief Sets/Unsets the internal memory data to be updated in the next assemble/contribute call */
    void SetUpdateMem(bool update = 1) override;
    
    /** @brief Gets the internal memory data to be updated in the next assemble/contribute call */
    bool GetUpdateMem() override;
    /** @}*/
protected:

    /** @brief Shared pointer to material memory */
    std::shared_ptr<TPZAdmChunkVector<TMem>> fMemory;

    /** @brief Default memory settings */
    TMem fDefaultMem;

    /** @brief Flag to indicate whether the memory data are to be updated in an assemble loop */
    bool fUpdateMem{false};
};

template<class TMem>
class TPZMatWithMemBC : public TPZMatWithMem<TMem>{
protected:
    // this method is your chance to verify if the material to which this
    // BC interface applies is compatible with this boundary interface
    // it is called in the method SetMaterial of class TPZBndCondBase
    void SetMaterialImpl(TPZMaterial* mat){}
    //do nothing
};

template <class TMem>
TPZMatWithMem<TMem>::TPZMatWithMem() :
fMemory(new TPZAdmChunkVector<TMem>())
{}

template <class TMem>
TPZMatWithMem<TMem>::TPZMatWithMem(const TPZMatWithMem<TMem> &mat) :
TPZRegisterClassId(&TPZMatWithMem::ClassId),
fMemory(new TPZAdmChunkVector<TMem>()), fDefaultMem(mat.fDefaultMem), fUpdateMem(mat.fUpdateMem) {
    *fMemory = *mat.fMemory;
}

template <class TMem>
TPZMatWithMem<TMem>&
TPZMatWithMem<TMem>::operator=(const TPZMatWithMem<TMem> &mat)
{
    fDefaultMem = mat.fDefaultMem;
    fUpdateMem = mat.fUpdateMem;
    if(!fMemory){
        fMemory(new TPZAdmChunkVector<TMem>());
    }
    *fMemory = *mat.fMemory;
}


template <class TMem>
void TPZMatWithMem<TMem>::Print(std::ostream &out) const {
    out << __PRETTY_FUNCTION__ << std::endl;
    out << "\nfDefaultMem = \n" << fDefaultMem;
    out << "\nfUpdateMem = " << fUpdateMem;
    int size = fMemory->NElements();
    out << "\nfMemory with " << size << " elements";
    for (int i = 0; i < size; i++) {
        out << "\nfMemory element : " << i << std::endl;
        this->MemItem(i).Print(out);
    }

}

template <class TMem>
void TPZMatWithMem<TMem>::PrintMem(std::ostream &out, const int memory) const {
    int size = fMemory->NElements();
    if (memory >= 0 && memory < size) {
        out << "fMemory element : " << memory << std::endl;
        this->MemItem(memory).Print(out);
    } else {
        out << "Memory index out of range : memory " << memory << " no elements " << size << std::endl;
    }
}

template <class TMem>
TMem &TPZMatWithMem<TMem>::MemItem(const int i) const {
    return fMemory.get()->operator [](i);
}

template <class TMem>
int TPZMatWithMem<TMem>::ClassId() const {
    return Hash("TPZMatWithMem") ^
        TMem().ClassId() << 1;
}

template <class TMem>
void TPZMatWithMem<TMem>::Write(TPZStream &buf, int withclassid) const {
    int updatemem = fUpdateMem;
    buf.Write(&updatemem);
    fDefaultMem.Write(buf, 0);
    TPZPersistenceManager::WritePointer(fMemory.get(), &buf);
}

template <class TMem>
void TPZMatWithMem<TMem>::Read(TPZStream &buf, void *context) {
    int updatemem;
    buf.Read(&updatemem);
    if (updatemem) {
        fUpdateMem = true;
    } else {
        fUpdateMem = false;
    }
    fDefaultMem.Read(buf, 0);
    fMemory = std::dynamic_pointer_cast<TPZAdmChunkVector<TMem> >(TPZPersistenceManager::GetSharedPointer(&buf));

}

template <class TMem>
std::shared_ptr<TPZAdmChunkVector<TMem>> & TPZMatWithMem<TMem>::GetMemory() {
    return fMemory;
}

template <class TMem>
void TPZMatWithMem<TMem>::SetMemory(std::shared_ptr<TPZAdmChunkVector<TMem>> & memory) {
    fMemory = memory;
}

template <class TMem>
int TPZMatWithMem<TMem>::PushMemItem(int sourceIndex) {
    int index = fMemory->AllocateNewElement();
    if (sourceIndex < 0) {
        this->ResetMemItem(index);
    } else {
        this->MemItem(index) = this->MemItem(sourceIndex);
    }
    return index;
}

template <class TMem>
void TPZMatWithMem<TMem>::FreeMemItem(int index) {
    fMemory->SetFree(index);
}

template <class TMem>
void TPZMatWithMem<TMem>::SetDefaultMem(TMem & defaultMem) {
    fDefaultMem = defaultMem;
}

template <class TMem>
void TPZMatWithMem<TMem>::SetUpdateMem(bool update) {
    fUpdateMem = update;
}

template <class TMem>
bool TPZMatWithMem<TMem>::GetUpdateMem() {
    return fUpdateMem;
}

#endif
