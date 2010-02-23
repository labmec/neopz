#ifndef PZCHECKCONSISTENCY
#define PZCHECKCONSISTENCY
/*
 *  pzcheckconsistency.h
 *  IP3D_v3
 *
 *  Created by Philippe Devloo on 27/12/09.
 *  Copyright 2009 UNICAMP. All rights reserved.
 *
 */
#include "pzsave.h"

#include <string>

/// class which implements an interface to check the consistency of two implementations
class TPZCheckConsistency
{
	/// a counter which will be used to compose the file name
	int fCounter;
	/// path where the file will be stored. If CHECKPATH is defined then this will be used as the file path, else the current directory will be the path
	std::string fPath;
	/// base file name. This name has to be unique to avoid overwrites
	std::string fFileName;
	/// boolean indicating whether the objects should be overwritten or not
	/**
	 * this flag is set to false by default
	 */
	bool fOverWrite;
	
	/// boolean indicating whether the object will write the objects to disk or read them
	/**
	 * the value of this variable is set default to true or false depending on the compiler directive READCHECK or WRITECHECK
	 */
	bool fWriteFlag;
	
public:
	/// constructor indicating the root filename
	TPZCheckConsistency(const std::string &filename);
	
	/// set the overwrite flag
	void SetOverWrite(bool flag = true);
	
	/// overwrite the read or write mode
	void SetReadMode()
	{
		fWriteFlag = false;
	}
	/// overwrite the read or write mode
	void SetWriteMode()
	{
		fWriteFlag = true;
	}
	/// this method goes reads or writes depending on the compiler directive
	/**
	 * when writing the method will write a binary copy of the object to a binary file
	 * when reading the method will read an object from the binary file and compare both copies
	 */
	bool CheckObject(TPZSaveable &obj);
};


#endif