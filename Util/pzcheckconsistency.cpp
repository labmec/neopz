/**
 * @file
 * @brief Contains the implementation of the TPZCheckConsistency methods. 
 * @author Created by Philippe Devloo on 27/12/09.
 */
/**
 * Copyright 2009 UNICAMP. All rights reserved.
 */

#include "pzcheckconsistency.h"
#include "pzbfilestream.h"
#include <sstream>
#include <fstream>

/// constructor indicating the root filename
TPZCheckConsistency::TPZCheckConsistency(const std::string &filename)
{
	fCounter = 0;
#ifdef CHECKPATH
	fPath = CHECKPATH;
#else
	fPath = "./";
#endif
	fOverWrite = false;
	fFileName = filename;
	fWriteFlag = true;
#ifdef READCHECK
	fWriteFlag = false;
#endif
#ifdef WRITECHECK
	fWriteFlag = true;
#endif
}

/// set the overwrite flag
void TPZCheckConsistency::SetOverWrite(bool flag)
{
	fOverWrite = flag;
}
/// this method goes reads or writes depending on the compiler directive
/**
 * when writing the method will write a binary copy of the object to a binary file
 * when reading the method will read an object from the binary file and compare both copies
 */
bool TPZCheckConsistency::CheckObject(TPZSaveable &obj)
{
	std::stringstream sout;
	sout << fCounter++;
	std::string name = fPath + fFileName + sout.str() + ".chk";
	TPZBFileStream file;
	if(!fWriteFlag)
	{
		file.OpenRead(name);
		TPZSaveable *copy = TPZSaveable::Restore(file,NULL);
		bool result = obj.Compare(copy,fOverWrite);
		delete copy;
		return result;
	}
	else
	{
		file.OpenWrite(name);
		obj.Write(file,true);
		return true;
	}
	//	std::ofstream file(name.str(),std::ios_base::binary);
}
