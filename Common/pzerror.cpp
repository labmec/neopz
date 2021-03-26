#include "pzerror.h"

void DebugStopImpl(const char *fileName, const std::size_t lineN)
{
#ifdef WIN32
	//ShowMessage("Erro encontrado! Entre em contato com o suporte do programa!");
#endif
#ifdef PZDEBUG
	std::cerr << "\n\nYour chance to put a breakpoint at " << fileName<< ":"<< lineN <<  "\n";
#endif
                                                                            throw std::bad_exception();
	
}
