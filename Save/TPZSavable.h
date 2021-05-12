/**
 * @file
 * @brief Contains declaration of the TPZSavable class which defines the interface to save and restore objects from TPZStream objects.
 */

#ifndef TPZSAVABLE_H
#define TPZSAVABLE_H

#include <string>     // for string
#include <set>        // for std::set
#include <map>        // for std::map
#include <utility>    // for pair
#include <list>       // for std::list
#include "pzerror.h"  // for DebugStop
class TPZSavable;
class TPZStream;
class TPZRestoreClassBase;
class TPZChunkTranslator;
/**
 * \addtogroup save
 * @{
 */
/** @brief Identifier as saveable object */
const int TPZSAVEABLEID = -1;

/** @brief Typedef of TPZRestore_t */
typedef TPZRestoreClassBase * TPZRestore_t;


/** 
 * The SAVEABLE NOTES are used to debug read and write operations on
 * TPZSavable objects. Do not use it for other purposes.
 * To enable it, define the DEBUG_SAVEABLE macro.
 */
//#define DEBUG_SAVEABLE

#ifdef PZDEBUG_SAVEABLE
#define SAVEABLE_STR_NOTE(buf,str) { std::string msg(str); buf.Write(&msg,1); }
#define SAVEABLE_SKIP_NOTE(buf) { std::string str; buf.Read(&str,1); }
#else
#define SAVEABLE_STR_NOTE(buf,str)
#define SAVEABLE_SKIP_NOTE(buf)
#endif


//Change the default constructor of TPZRegisterClassId to
// TPZRegisterClassId() = delete;
// in order to check all TPZSavable children classes for ClassId()
class TPZRegisterClassId {
public:
    // this matches the signature of 'ClassId() const'
    template <typename T>
    TPZRegisterClassId(int (T::*)() const) {
    }
    TPZRegisterClassId() = default;
};

/**
 * @brief This class defines the interface to save and restore objects from TPZStream objects. \ref save "Persistency"
 */
/**
 * This class defines the interface a class needs to implement (override) in order to become persistent
 * Several static utility methods have been defined to make saving and restoring of vectors of
 * objects easier
 */
class TPZSavable : public virtual TPZRegisterClassId {
	
public :
	/** @brief This static function guarantees that the gMap object is available when needed */
	static std::set<TPZRestoreClassBase*> &RestoreClassSet() {
		static std::set<TPZRestoreClassBase*> gRestoreClassSet;
		return gRestoreClassSet;
	}
        
	/** @brief This static function guarantees that the gMap object is available when needed */
	static std::map<int,TPZRestore_t> &ClassIdMap() {
		static std::map<int,TPZRestore_t> gClassIdMap;
		return gClassIdMap;
	}
	
	
public:
	
    TPZSavable() : TPZRegisterClassId(&TPZSavable::ClassId) { } 
	virtual ~TPZSavable()
	{
	}
	
	/** @brief Define the class id associated with the class */
	/**
	 * This id has to be unique for all classes
	 * A non unique id is flagged at the startup of the program
	 */
	public:
        
	virtual std::list<std::map<std::string, uint64_t>> VersionHistory() const;
	virtual std::pair<std::string, uint64_t> Version() const;
    
    static std::pair<std::string, uint64_t> NeoPZVersion();

    
    /** @name Read and Write methods
     *  Methods used for reading and writing to file.
     */
    //@{
    /** @brief Unique identifier for each class */
    virtual int ClassId() const = 0;
	/** @brief Writes this object to the TPZStream buffer. Include the classid if `withclassid = true` */
	virtual void Write(TPZStream &buf, int withclassid) const;
	
	/** @brief Reads an objects from the TPZStream buffer. */
	virtual void Read(TPZStream &buf, void *context);
    //@}
	/** @brief Compares the object for identity with the object pointed to, eventually copy the object */
	/**
	 * Compares both objects bitwise for identity. Put an entry in the log file if different
	 * overwrite \n the calling object if the override flag is true
	 */
	virtual bool Compare(TPZSavable *copy, bool override = false);
	
	/** @brief Compares the object for identity with the object pointed to, eventually copy the object */
	/**
	 * Compares both objects bitwise for identity. Put an entry in the log file if different
	 * generate \n an interrupt if the override flag is true
	 */
	virtual bool Compare(TPZSavable *copy, bool override = false) const;

        static void Register(TPZRestoreClassBase *restore);
	
        static void RegisterClassId(int classid, TPZRestore_t fun);
	
        static TPZSavable *CreateInstance(const int &classId);
};


class TPZRestoreClassBase{
public :
    virtual TPZSavable *Restore()=0;
    virtual TPZChunkTranslator *GetTranslator()=0;
};

/**
 * @brief Implements an interface to register a class id and a restore function. \ref save "Persistence"
 */

/**
 * A declaration of the type "template class<classname> put in .cpp file does the trick \n
 * The static object which is "automatically" created calls the proper interface of the TPZSavable class
 */
template<class T>
class TPZRestoreClass : public TPZRestoreClassBase {
public:

    /** @brief Constructor */
    TPZRestoreClass() {
//        std::cout << __PRETTY_FUNCTION__ << "Registering\n";
        TPZSavable::Register(this);
    }

    /** @brief Restores object from Map based in classid into the buf */
    virtual TPZSavable *Restore() override {
        T *ptr = new T;
        return ptr;
    }
    
    virtual TPZChunkTranslator *GetTranslator() override {
        return nullptr;
    }
    
private:
    static TPZRestoreClass<T> gRestoreObject;
};

template<class T>
TPZRestoreClass<T> TPZRestoreClass<T>::gRestoreObject;


template<class T, class TranslatorType>
class TPZRestoreClassWithTranslator : public TPZRestoreClassBase {
public:

    /** @brief Constructor */
    TPZRestoreClassWithTranslator() {
        TPZSavable::Register(this);
    }
    
    ~TPZRestoreClassWithTranslator();

    /** @brief Restores object from Map based in classid into the buf */
    virtual TPZSavable *Restore() override {
        T *ptr = new T;
        return ptr;
    }
    
    virtual TPZChunkTranslator *GetTranslator() override {
        if (!gTranslator){
            gTranslator = new TranslatorType();
        }
        return gTranslator;
    }
    
private:
    static TPZRestoreClassWithTranslator<T,TranslatorType> gRestoreObject;
    static TPZChunkTranslator *gTranslator;
};

#include "TPZChunkTranslator.h"
template<class T, class TranslatorType>
TPZRestoreClassWithTranslator<T, TranslatorType>::~TPZRestoreClassWithTranslator() {
    if (gTranslator) {
        delete gTranslator;
    }
}

template<class T, class TranslatorType>
TPZRestoreClassWithTranslator<T,TranslatorType> TPZRestoreClassWithTranslator<T,TranslatorType>::gRestoreObject;

template<class T, class TranslatorType>
TPZChunkTranslator *TPZRestoreClassWithTranslator<T,TranslatorType>::gTranslator = NULL;

///@}
#endif //TPZSAVABLE_H
