/**
 * @file
 * @brief Defines functions to help with invalid arithmetic operations.
 *
 *   This file implements functions previously implemented in http://www-personal.umich.edu/~williams/archive/computation/fe-handling-example.c in order to improve portability of NeoPZ code and allow feenableexcept function to be called in a macOS environment. These functions were slightly modified.
 */
#ifndef __FPO_EXCEPTIONS_H
#define __FPO_EXCEPTIONS_H

#ifdef WIN32

#include <float.h>

#ifndef _EM_OVERFLOW
#define _EM_OVERFLOW EM_OVERFLOW
#endif//_EM_OVERFLOW

#ifndef _EM_UNDERFLOW
#define _EM_UNDERFLOW EM_UNDERFLOW
#endif//_EM_UNDERFLOW

#ifndef _EM_INVALID
#define _EM_INVALID EM_INVALID
#endif//_EM_INVALID

#ifndef _EM_ZERODIVIDE
#define _EM_ZERODIVIDE EM_ZERODIVIDE
#endif//_EM_ZERODIVIDE

#else 

#include <fenv.h>

#ifdef MACOSX


inline int
fegetexcept (void)
{
    static fenv_t fenv;
    
    return fegetenv (&fenv) ? -1 : (fenv.__control & FE_ALL_EXCEPT);
}

inline int
feenableexcept (unsigned int excepts)
{
    static fenv_t fenv;
    unsigned int new_excepts = excepts & FE_ALL_EXCEPT,
    old_excepts;  // previous masks
    
    if ( fegetenv (&fenv) ) return -1;
    old_excepts = fenv.__control & FE_ALL_EXCEPT;
    
    // unmask
    fenv.__control &= ~new_excepts;
    fenv.__mxcsr   &= ~(new_excepts << 7);
    
    return ( fesetenv (&fenv) ? -1 : old_excepts );
}

inline int
fedisableexcept (unsigned int excepts)
{
    static fenv_t fenv;
    unsigned int new_excepts = excepts & FE_ALL_EXCEPT,
    old_excepts;  // all previous masks
    
    if ( fegetenv (&fenv) ) return -1;
    old_excepts = fenv.__control & FE_ALL_EXCEPT;
    
    // mask
    fenv.__control |= new_excepts;
    fenv.__mxcsr   |= new_excepts << 7;
    
    return ( fesetenv (&fenv) ? -1 : old_excepts );
}


// x87 fpu
#define getx87cr(x)    asm ("fnstcw %0" : "=m" (x));
#define setx87cr(x)    asm ("fldcw %0"  : "=m" (x));
#define getx87sr(x)    asm ("fnstsw %0" : "=m" (x));

// SIMD, gcc with Intel Core 2 Duo uses SSE2(4)
#define getmxcsr(x)    asm ("stmxcsr %0" : "=m" (x));
#define setmxcsr(x)    asm ("ldmxcsr %0" : "=m" (x));

#include <signal.h>
#include <stdio.h>   // printf()
#include <stdlib.h>  // abort(), exit()

static std::string fe_code_name[] = {
    "FPE_NOOP",
    "FPE_FLTDIV", "FPE_FLTINV", "FPE_FLTOVF", "FPE_FLTUND",
    "FPE_FLTRES", "FPE_FLTSUB", "FPE_INTDIV", "FPE_INTOVF"
    "FPE_UNKNOWN"
};

/* SAMPLE ALTERNATE FP EXCEPTION HANDLER
 
 The sample handler just reports information about the
 exception that invoked it, and aborts.  It makes no attempt
 to restore state and return to the application.
 
 More sophisticated handling would have to confront at least
 these issues:
 
 * interface to the system context for restoring state
 * imprecision of interrupts from hardware for the intel x87
 fpu (but not the SIMD unit, nor the ppc)
 * imprecision of interrupts from system software
 */
inline void
fhdl ( int sig, siginfo_t *sip, ucontext_t *scp )
{
    int fe_code = sip->si_code;
    unsigned int excepts = fetestexcept (FE_ALL_EXCEPT);
    
    switch (fe_code)
    {
#ifdef FPE_NOOP  // occurs in OS X
        case FPE_NOOP:   fe_code = 0; break;
#endif
        case FPE_FLTDIV: fe_code = 1; break; // divideByZero
        case FPE_FLTINV: fe_code = 2; break; // invalid
        case FPE_FLTOVF: fe_code = 3; break; // overflow
        case FPE_FLTUND: fe_code = 4; break; // underflow
        case FPE_FLTRES: fe_code = 5; break; // inexact
        case FPE_FLTSUB: fe_code = 6; break; // invalid
        case FPE_INTDIV: fe_code = 7; break; // overflow
        case FPE_INTOVF: fe_code = 8; break; // underflow
        default: fe_code = 9;
    }
    
    if ( sig == SIGFPE )
    {
        unsigned short x87cr,x87sr;
        unsigned int mxcsr;
        
        getx87cr (x87cr);
        getx87sr (x87sr);
        getmxcsr (mxcsr);
        printf ("X87CR:   0x%04X\n", x87cr);
        printf ("X87SR:   0x%04X\n", x87sr);
        printf ("MXCSR:   0x%08X\n", mxcsr);

        
        printf ("signal:  SIGFPE with code %s\n", fe_code_name[fe_code].c_str());
        printf ("invalid flag:    0x%04X\n", excepts & FE_INVALID);
        printf ("divByZero flag:  0x%04X\n", excepts & FE_DIVBYZERO);
    }
    else printf ("Signal is not SIGFPE, it's %i.\n", sig);
    abort();
}
#endif //MACOSX
//#endif //DEBUG

#endif //WIN32


struct TExceptionManager {
	private:
#ifdef WIN32
	unsigned int fPrevConfig;
#else
    fenv_t fPrevConfig;
#endif //WIN32
public:
    TExceptionManager(){
#ifdef WIN32
        _controlfp_s(&fPrevConfig, 0, 0);//saves current state of fpu
        _controlfp(1, _EM_OVERFLOW);
        _controlfp(1, _EM_INVALID);
        _controlfp(1, _EM_ZERODIVIDE);
#else
        fegetenv(&fPrevConfig);//saves current state of fpu
        feraiseexcept(FE_INVALID | FE_OVERFLOW | FE_DIVBYZERO);
#endif //WIN32
    }
    
    ~TExceptionManager(){
#ifdef WIN32
        unsigned int temp;
        _controlfp_s(&temp, fPrevConfig, _MCW_EM);//restores previous sates of fpu
#else
        fesetenv(&fPrevConfig);//restores previous sates of fpu
#endif //WIN32
    }
};


#endif //__FPO_EXCEPTIONS_H
