#include "pzerror.h"
#include <iostream>

#include <stdexcept>
void pzinternal::DebugStopImpl(const char *fileName, const std::size_t lineN)
{
#ifdef WIN32
	//ShowMessage("Erro encontrado! Entre em contato com o suporte do programa!");
#endif
    //#ifdef PZDEBUG
	PZError << "Your chance to put a breakpoint at " << fileName<< ":"<< lineN <<  "\n";
    //#endif
                                                                            throw std::bad_exception();
	
}
