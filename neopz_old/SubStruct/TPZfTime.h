#ifndef TPZFTIME_H
#define TPZFTIME_H

/*
 *  ftime.h
 *  SubStruct
 *
 *  Created by Bandit on 7/20/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */


#include <iostream>
#include <sys/timeb.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <sstream>

/// Calculate the Times
class TPZfTime
{
	
public:

	/// Start the timer when the object is created
	TPZfTime();
	
	~TPZfTime();
	
	/// When called, returns the time since the creation of the object in a string
	std::string ReturnTimeString();
	
	/// When called, returns the time since the creation of the object in a double
	double ReturnTimeDouble();
	
	
	
	
private:

	struct timeb finicio, ffinal;

	
};

#endif