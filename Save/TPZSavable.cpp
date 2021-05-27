/**
 * @file
 * @brief Contains the TPZSavable methods.
 */

#include "TPZSavable.h"
#include <iostream>                        // for operator<<, basic_ostream
#include <vector>                          // for allocator
#include "TPZStream.h"                     // for TPZStream
#include "pzlog.h"                         // for glogmutex, LOGPZ_ERROR
#include <list>

#include "pzlog.h"
#include "TPZPersistenceManager.h"
#include "TPZChunkTranslator.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.savable");
static TPZLogger loggerCheck("pz.checkconsistency");
#endif

std::list<std::map<std::string, uint64_t>> TPZSavable::VersionHistory() const {
    std::list<std::map<std::string, uint64_t>> history;
    std::map<std::string, uint64_t> versionMap;
    versionMap.insert(std::make_pair("NeoPZ", 1));
    history.push_back(versionMap);
    versionMap.clear();
    versionMap.insert(std::make_pair("NeoPZ", 2));
    history.push_back(versionMap);
    versionMap.clear();
    versionMap.insert(std::make_pair("NeoPZ", 3));
    history.push_back(versionMap);
    versionMap.clear();
    versionMap.insert(std::make_pair("NeoPZ", 4));
    history.push_back(versionMap);
    return history;
}

std::pair<std::string, uint64_t> TPZSavable::Version() const {
    return TPZSavable::NeoPZVersion();
}

std::pair<std::string, uint64_t> TPZSavable::NeoPZVersion() {
    return std::make_pair("NeoPZ", 4);
}

void TPZSavable::Write(TPZStream &buf, int withclassid) const
{
    DebugStop();
	if(withclassid) { 
		int var = ClassId();
		if(var == -1)
		{
			std::cout << "TPZSavable::Write const with classid -1 expect trouble\n";
		}
		//buf.Write(&var,1);
	}
}


void TPZSavable::Read(TPZStream &buf, void *context)
{}

void TPZSavable::Register(TPZRestoreClassBase *restore) {
	std::set<TPZRestoreClassBase*>::iterator it;
	it = RestoreClassSet().find(restore);
	if(it != RestoreClassSet().end()) 
	{
		std::cout << "TPZSavable::Register duplicate RestoreClass " << std::endl;
                DebugStop();
	}
	RestoreClassSet().insert(restore);
#ifdef PZ_LOG
//    std::cout << "Registering object " << typeid(restore).name() <<
//        std::endl;
#endif
}


void TPZSavable::RegisterClassId(int classid, TPZRestore_t fun) 
{

  std::map<int,TPZRestore_t>::iterator it;
	it = ClassIdMap().find(classid);
	if(it != ClassIdMap().end()) 
	{
    auto itsrest = it->second->Restore();
    auto funrest = fun->Restore();
		std::cout << "TPZSavable::Register duplicate classid\n"
              << "class id : "<<itsrest->ClassId()<<'\n'
              << "type id : " << typeid(*(itsrest)).name() <<'\n'
              << "conflicted with:\n"
              << "class id : "<<classid<<'\n'
              << "type id : " << typeid(*(funrest)).name() <<std::endl;
    DebugStop();
	}
	ClassIdMap()[classid] = fun;
    if (fun->GetTranslator()){
        fun->GetTranslator()->SetClassId(classid);
    }

}

/// Compare the object for identity with the object pointed to, eventually copy the object
/**
 * compare both objects bitwise for identity. Put an entry in the log file if different
 * overwrite the calling object if the override flag is true
 */
bool TPZSavable::Compare(TPZSavable *copy, bool override)
{
	std::stringstream sout;
	sout << "Class id " << ClassId() << " Compare needs to be implemented";
	LOGPZ_ERROR(loggerCheck,sout.str())
	return false;
}

/// Compare the object for identity with the object pointed to, eventually copy the object
/**
 * compare both objects bitwise for identity. Put an entry in the log file if different
 * overwrite the calling object if the override flag is true
 */
bool TPZSavable::Compare(TPZSavable *copy, bool override) const
{
	std::stringstream sout;
	sout << "Class id " << ClassId() << " Compare needs to be implemented";
	LOGPZ_ERROR(loggerCheck,sout.str())
	return false;
}


TPZSavable *TPZSavable::CreateInstance(const int &classId) {
    
    std::map<int,TPZRestore_t>::const_iterator it;
    it = ClassIdMap().find(classId);
    if(it == ClassIdMap().end()) {
        std::cout << "TPZSavable trying to restore unknown object with classId " << classId << std::endl;
#ifdef PZ_LOG
        {
            std::stringstream sout;
            sout << __PRETTY_FUNCTION__ << " trying to restore unknown object with classId " << classId;
            LOGPZ_ERROR(logger,sout.str().c_str());
        }
#endif
        TPZPersistenceManager::CloseRead();
        DebugStop();
    }
    
    TPZRestore_t fun= it->second;
    return fun->Restore();
}
