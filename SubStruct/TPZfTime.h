#include <iostream>
#include <sys/timeb.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <strstream>

class TPZfTime
{

public:

TPZfTime();

~TPZfTime();

std::strstream ReturnTimeString();
double ReturnTimeDouble();


};
