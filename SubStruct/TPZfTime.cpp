/*
 *  ftime.cpp
 *  SubStruct
 *
 *  Created by Bandit on 7/20/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "TPZftime.h"
#include <iostream>
#include <sys/timeb.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <sstream>

using namespace std;

TPZfTime::TPZfTime()
{
	ftime(&finicio);
}

TPZfTime::~TPZfTime()
{
}

string TPZfTime::ReturnTimeString()
{ 

	ftime(&ffinal);
	double time = ((double) ffinal.time + ((double) ffinal.millitm * 0.001)) - ((double) finicio.time + ((double) finicio.millitm * 0.001));
	stringstream oss;
	string str;
	
	oss << time << " seconds\n";
	str = oss.str();
	
   // printf("%f seconds\n", ((double) ffinal.time + ((double) ffinal.millitm * 0.001)) - ((double) finicio.time + ((double) finicio.millitm * 0.001)));
	return str;
}

double TPZfTime::ReturnTimeDouble()
{
	ftime(&ffinal);
	double time = ((double) ffinal.time + ((double) ffinal.millitm * 0.001)) - ((double) finicio.time + ((double) finicio.millitm * 0.001));
	
	return time;
}
