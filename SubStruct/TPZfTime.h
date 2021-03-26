/**
 * @file
 * @brief Contains the TPZfTime class which calculates times.
 * @since 7/20/10.
 */

#ifndef TPZFTIME_H
#define TPZFTIME_H

#include <iostream>
#include <sys/timeb.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <sstream>

/** 
 * @brief Calculate the Times. \ref util "Utility"
 * @ingroup util
 */
class TPZfTime
{
	
public:
	
	/** @brief Start the timer when the object is created */
	TPZfTime();
	/** @brief Default destructor */
	~TPZfTime();
	
	/** @brief When called, returns the time since the creation of the object in a string */
	std::string ReturnTimeString();
	
	/** @brief When called, returns the time since the creation of the object in a double */
	double ReturnTimeDouble();
	
private:
	/** @brief Initial and final time calculates. */
	struct timeb finicio, ffinal;
	
};

#endif
