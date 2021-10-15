#include "pzerror.h"
#include <iostream>

#include <stdexcept>

#include <execinfo.h>
#include <stdio.h>
void pzinternal::DebugStopImpl(const char *fileName, const std::size_t lineN)
{
#ifdef WIN32
	//ShowMessage("Erro encontrado! Entre em contato com o suporte do programa!");
#endif
    //#ifdef PZDEBUG
	PZError << "\n\nYour chance to put a breakpoint at " << fileName<< ":"<< lineN <<  "\n";
    //#endif

    void* callstack[128];
    int i, frames = backtrace(callstack, 128);
    char** strs = backtrace_symbols(callstack, frames);
    for (i = 0; i < frames; ++i) {
        printf("%s\n", strs[i]);
    }
    free(strs);


    throw std::bad_exception();
	
}
