/**
 * @file
 * @brief Contains declaration of the TPZSaveable class which defines the interface to save and restore objects from TPZStream objects.
 */

#ifndef TPZSAVEABLE_H
#define TPZSAVEABLE_H

#include <map>
#include <vector>
#include <set>


#include "pzvec.h"
#include "pzmanvector.h"
#include "pzadmchunk.h"
#include "pzreal.h"
#include "tpzautopointer.h"

class TPZSaveable;
class TPZStream;
/**
 * \addtogroup save
 * @{
 */
/** @brief Identifier as saveable object */
const int TPZSAVEABLEID = -1;

/** @brief Typedef of TPZRestore_t */
typedef TPZSaveable *(*TPZRestore_t)(TPZStream &,void *);


/** 
 * The SAVEABLE NOTES are used to debug read and write operations on
 * TPZSaveable objects. Do not use it for other purposes.
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

/**
 * @brief This class defines the interface to save and restore objects from TPZStream objects. \ref save "Persistency"
 */
/**
 * This class defines the interface a class needs to implement (override) in order to become persistent
 * Several static utility methods have been defined to make saving and restoring of vectors of
 * objects easier
 */
class TPZSaveable {
	
#ifndef ELLIPS
	
	/** @brief This static function garantes that the gMap object is available when needed */
	static std::map<int,TPZRestore_t> &Map() {
		static std::map<int,TPZRestore_t> gMap;
		return gMap;
	}
#endif
	
public:
    
#ifndef ELLIPS
    
    static std::map<std::string, const int> &versionInfo() {
        static std::map<std::string, const int> vInfoMap;
        vInfoMap.insert(std::pair<std::string, const int>("NeoPZ", 1));//will be inserted only once
        return vInfoMap;
    }
    
#endif
	
	virtual ~TPZSaveable()
	{
	}
	
	/** @brief Define the class id associated with the class */
	/**
	 * This id has to be unique for all classes
	 * A non unique id is flagged at the startup of the program
	 */
	virtual int ClassId() const ;
	
	/** @brief Writes this object to the TPZStream buffer. Include the classid if withclassid = true */
	virtual void Write(TPZStream &buf, int withclassid);
	
	/** @brief Writes this object to the TPZStream buffer. Include the classid if withclassid = true */
	virtual void Write(TPZStream &buf, int withclassid) const;
	
	/** @brief read objects from the stream */
	virtual void Read(TPZStream &buf, void *context);
	
	/** @brief Compares the object for identity with the object pointed to, eventually copy the object */
	/**
	 * Compares both objects bitwise for identity. Put an entry in the log file if different
	 * overwrite \n the calling object if the override flag is true
	 */
	virtual bool Compare(TPZSaveable *copy, bool override = false);
	
	/** @brief Compares the object for identity with the object pointed to, eventually copy the object */
	/**
	 * Compares both objects bitwise for identity. Put an entry in the log file if different
	 * generate \n an interrupt if the override flag is true
	 */
	virtual bool Compare(TPZSaveable *copy, bool override = false) const;
	
	static void Register(int classid, TPZRestore_t fun);
	
	static TPZSaveable *Restore(TPZStream &buf, void *context);
    
    virtual void AssignPointers(const TPZVec<int> & idsVec){
        DebugStop();
    }
    
};

/** @brief Restores object from Map, classid is in buf */
template<class T>
TPZSaveable *Restore(TPZStream &buf, void *context) {
	T *ptr = new T;
	ptr->Read(buf,context);
	return ptr;
}

#ifndef ELLIPS

/**
 * @brief Implements an interface to register a class id and a restore function. \ref save "Persistency"
 */
/**
 * A declaration of the type "template class<classname, classid> put in .cpp file does the trick \n
 * The static object which is "automatically" created calls the proper interface of the TPZSaveable class
 */
template<class T, int N>
class TPZRestoreClass {
public:
	/** @brief Constructor */
	TPZRestoreClass()
	{
#ifdef PZDEBUG 
		std::string func_name = __PRETTY_FUNCTION__;
#ifndef WIN32
		std::cout << func_name << std::endl;
#endif
#endif
		TPZSaveable::Register(N,Restore);
	}
public:
	/** @brief Restores object from Map based in classid into the buf */
	static TPZSaveable *Restore(TPZStream &buf, void *context) {
		T *ptr = new T;
		ptr->Read(buf,context);
		return ptr;
	}
private:
	
	static TPZRestoreClass gRestoreObject;
};

template<>
inline TPZSaveable *TPZRestoreClass<TPZSaveable,-1>::Restore(TPZStream &buf, void *context)
{
	return 0;
}
template<class T, int N>
TPZRestoreClass<T,N> TPZRestoreClass<T,N>::gRestoreObject;

/// To restore object
template<>
inline TPZSaveable *Restore<TPZSaveable>(TPZStream &buf, void *context) {
	return 0;
}

/** @} */

#endif //ellips

//#include "TPZStream.h"
//
//
//TPZSaveable *TPZSaveable::Restore(TPZStream &buf, void *context) {
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
//        std::cout << "TPZSaveable trying to restore unknown object " << classid << std::endl;
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

#endif //TPZSAVEABLE_H
