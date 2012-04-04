/**
 * @file
 * @brief Contains the TPZSaveable methods.
 */

#include "pzsave.h"
#include "pzfilebuffer.h"

#include <iostream>
#include <fstream>
#include <sstream>


#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.saveable"));
static LoggerPtr loggerCheck(Logger::getLogger("pz.checkconsistency"));
#endif

using namespace std;


int TPZSaveable::ClassId() const {
	return -1;
}

void TPZSaveable::Write(TPZStream &buf, int withclassid) 
{
	if(withclassid) { 
		int var = ClassId();
		if(var == -1)
		{
			cout << "TPZSaveable::Write with classid -1 expect trouble\n";
		}
		buf.Write(&var,1);
	}
}

void TPZSaveable::Read(TPZStream &buf, void *context)
{}

void TPZSaveable::Register(int classid, TPZRestore_t fun) 
{
#ifndef ELLIPS
	map<int,TPZRestore_t>::iterator it;
	it = Map().find(classid);
	if(it != Map().end()) 
	{
		cout << "TPZSaveable::Register duplicate classid " << classid << endl;
		return;
	}
	Map()[classid] = fun;
#endif
}

TPZSaveable *TPZSaveable::Restore(TPZStream &buf, void *context) {
#ifndef ELLIPS
	int classid;
	buf.Read(&classid,1);
	map<int,TPZRestore_t>::iterator it;
	it = Map().find(classid);
	if(it == Map().end()) 
	{
		cout << "TPZSaveable trying to restore unknown object " << classid << endl;
		{
			std::stringstream sout;
			sout << __PRETTY_FUNCTION__ << " trying to restore unknown object " << classid;
			LOGPZ_ERROR(logger,sout.str().c_str());
		}
		return 0;
	}
	
	TPZRestore_t fun= it->second;
	return (*fun)(buf,context);
#else
	return 0;
#endif
}

/// Compare the object for identity with the object pointed to, eventually copy the object
/**
 * compare both objects bitwise for identity. Put an entry in the log file if different
 * overwrite the calling object if the override flag is true
 */
bool TPZSaveable::Compare(TPZSaveable *copy, bool override)
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
bool TPZSaveable::Compare(TPZSaveable *copy, bool override) const
{
	std::stringstream sout;
	sout << "Class id " << ClassId() << " Compare needs to be implemented";
	LOGPZ_ERROR(loggerCheck,sout.str())
	return false;
}

template class TPZRestoreClass<TPZSaveable, -1>;
