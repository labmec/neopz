/**
 * @file
 * @brief Contains declaration of the TPZSavable class which defines the interface to save and restore objects from TPZStream objects.
 */

#ifndef TPZSAVABLE_H
#define TPZSAVABLE_H

#include <string>     // for string
#include <set>        // for set
#include <map>        // for map
#include <utility>    // for pair
#include "pzerror.h"  // for DebugStop
#include "pzvec.h"    // for TPZVec
class TPZSavable;  // lines 20-20
class TPZStream;  // lines 21-21
class TPZRestoreClassBase;
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


//Search for CREATECLASSID to find classes with previously missing ClassId()
//Comment the default constructor of TPZRegisterClassId in order to 
//check all TPZSavable  children classes for ClassId()
class TPZRegisterClassId {
public:
    // this matches the signature of 'ClassId() const'
    template <typename T>
    TPZRegisterClassId(int (T::*)() const) {
    }
    TPZRegisterClassId() {}
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
	
#ifndef ELLIPS
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
	
#endif
	
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
virtual int ClassId() const = 0;

        
	virtual std::pair<std::string, uint64_t> Version() const;
        	
	/** @brief Writes this object to the TPZStream buffer. Include the classid if withclassid = true */
	//virtual void Write(TPZStream &buf, int withclassid) const;
	
	/** @brief Writes this object to the TPZStream buffer. Include the classid if withclassid = true */
	virtual void Write(TPZStream &buf, int withclassid) const;
	
	/** @brief read objects from the stream */
	virtual void Read(TPZStream &buf, void *context);
	
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

#ifndef ELLIPS

class TPZRestoreClassBase{
public :
    virtual TPZSavable *Restore()=0;
};

/**
 * @brief Implements an interface to register a class id and a restore function. \ref save "Persistence"
 */
/**
 * A declaration of the type "template class<classname, classid> put in .cpp file does the trick \n
 * The static object which is "automatically" created calls the proper interface of the TPZSavable class
 */
template<class T>
class TPZRestoreClass : public TPZRestoreClassBase {
public:
	/** @brief Constructor */
	TPZRestoreClass() {
#ifdef PZDEBUG 
		std::string func_name = __PRETTY_FUNCTION__;
#endif
		TPZSavable::Register(this);
	}
public:
	/** @brief Restores object from Map based in classid into the buf */
	virtual TPZSavable *Restore() {
		T *ptr = new T;
		return ptr;
	}
private:
	
	static TPZRestoreClass gRestoreObject;
};

template<class T>
TPZRestoreClass<T> TPZRestoreClass<T>::gRestoreObject;


/** @brief Restores object from Map, classid is in buf */
//template<class T>
//TPZSavable *Restore(TPZStream &buf, const int &id) {
//    T *ptr = new T;
//    //TODO:: ADD TO VECTOR
//    void *context = NULL;
//    ptr->Read(buf,context);
//    return ptr;
//}

/// To restore object
//template<>
//inline TPZSavable *Restore<TPZSavable>(TPZStream &buf, const int &id) {
//	return 0;
//}

/** @} */

#endif //ellips

//#include "TPZStream.h"
//
//
//TPZSavable *TPZSavable::Restore(TPZStream &buf, void *context) {
//    //#ifdef LOG4CXX
//    //    LoggerPtr logger(Logger::getLogger("pz.saveable"));
//    //    LoggerPtr loggerCheck(Logger::getLogger("pz.checkconsistency"));
//    //#endif
//    
//#ifndef ELLIPS
//    int classid;
//    buf.Read(&classid,1);
//    map<int,TPZRestore_t>::iterator it;
//    it = Map().find(classid);
//    if(it == Map().end())
//    {
//        std::cout << "TPZSavable trying to restore unknown object " << classid << std::endl;
//        {
//            std::stringstream sout;
//            sout << __PRETTY_FUNCTION__ << " trying to restore unknown object " << classid;
//#ifdef LOG4CXX
//            //LOGPZ_ERROR(logger,sout.str().c_str());
//#endif
//        }
//        return 0;
//    }
//    
//    TPZRestore_t fun= it->second;
//    return (*fun)(buf,context);
//#else
//    return 0;
//#endif
//}

#endif //TPZSAVABLE_H
