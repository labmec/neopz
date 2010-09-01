#ifndef TPZDATABASE_H
#define TPZDATABASE_H

/*
 *  TPZDataBase.h
 *  SubStruct
 *
 *  Created by Bandit on 8/9/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "TPZTimeTemp.h"
#include "pzvec.h"

class TPZDataBase
{
	
public:
	/// Constructor
	TPZDataBase();
	
	///TPZVec with all the read informations
	TPZVec <TPZTimeTemp> fTimes;
	
	/// Read de File saving the information in TPZVec fTimes
	void Read(std::string &FileName);
	
	/// Return the number of lines in the file fFileName
	int NumberOfLines();
	
	///Set File Name in fFileName
	void SetFileName(std::string &FileName);
	
private:
	/// File Name
	std::string fFileName;	
	
};

#endif