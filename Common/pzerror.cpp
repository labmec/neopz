#include "pzerror.h"

void pzinternal::DebugStopImpl(const char *fileName, const std::size_t lineN)
{
#ifdef WIN32
	//ShowMessage("Erro encontrado! Entre em contato com o suporte do programa!");
#endif
#ifdef PZDEBUG
	std::cerr << "Your chance to put a breakpoint at " << fileName<< ":"<< lineN <<  "\n";
#endif
                                                                            throw std::bad_exception();
	
}
