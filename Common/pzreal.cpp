/**
 * @file
 * @brief Contains the TPZCounter methods and the DebugStop() function.
 */

#include "pzreal.h"
#include <string>
#include <exception>

using namespace std;

#ifndef ELLIPS

TPZCounter TPZFlopCounter::gCount;

static string names[] = {
	"sum",
	"prod",
	"div",
	"sqrt",
	"pow",
	"cos",
	"sin",
	"acos",
	"asin",
	"atan",
	"exp",
	"log"
};

void TPZCounter::Print(std::ostream &out) const
{
	int i;
	for(i=0; i<gNumOp; i++)
    {
		out << names[i] << " " << fCount[i] << endl;
    }
}
std::ostream &operator<<(std::ostream &out,const TPZCounter &count)
{
	int i;
	for(i=0; i<gNumOp; i++)
    {
		out << names[i] << "/" << count.fCount[i] << "\t";
    }
	return out;
}

#endif

#ifdef WIN32
//#include <Dialogs.hpp>
#endif

void DebugStop()
{
#ifdef WIN32
	//ShowMessage("Erro encontrado! Entre em contato com o suporte do programa!");
#endif
	std::cout << "Your chance to put a breakpoint at " << __FILE__ <<  "\n";
	std::bad_exception myex;
	throw myex;
	
}

/**
 * Function erf (Error function) implemented in 
 * http://www.johndcook.com/cpp_erf.html
 */
REAL erf(REAL arg) {
	REAL a1 = 0.254829592;
	REAL a2 = -0.284496736;
	REAL a3 = 1.421413741;
	REAL a4 = -1.453152027;
	REAL a5 = 1.061405429;
	REAL p = 0.3275911;
	// Save the sign of arg
	int sign = 1;
	if(arg < 0) sign = -1;
	arg = fabs(arg);

	// A&S formula 7.1.26
	REAL t = 1.0/(1.0 + p*arg);
	REAL y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-arg*arg);

	return sign*y;
}
