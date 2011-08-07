/**
 * @file
 * @brief Contains the TPZCounter methods and the DebugStop() function.
 */
//
// C++ Implementation: pzreal
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2004
//
// Copyright: See COPYING file that comes with this distribution
//

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
#include <Dialogs.hpp>
#endif

void DebugStop()
{
#ifdef WIN32
	ShowMessage("Erro encontrado! Entre em contato com o suporte do programa!");
#endif
	std::cout << "Your chance to put a breakpoint at " << __FILE__ <<  "\n";
	std::bad_exception myex;
	throw myex;
	
}
