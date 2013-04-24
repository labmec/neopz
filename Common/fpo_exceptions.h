/**
 * @file
 * @brief Defines functions to help with invalid arithmetic operations
 */
#ifndef __FPO_EXCEPTIONS_H
#define __FPO_EXCEPTIONS_H

#ifndef WIN32

#ifdef DEBUG
#include <iostream>
#include <stdlib.h>
#include <xmmintrin.h>
#include <signal.h>
#include <fenv.h>

/** 
 * @ingroup common
 * @brief Declares a handler function to deal with SIGFPE exception
 */
void InvalidFPOHandler(int signo) {
    switch(signo) {
        case SIGFPE: std::cout << "ERROR : Invalid Arithmetic operation." << std::endl; break;
    }
    exit(signo);
}

/** 
 * @ingroup common
 * @brief This macro enable exceptions during invalid FP operations. (Use this macro inside the main function)
 */
#define ENABLE_FPO_EXCEPTIONS _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);

/** 
 * @ingroup common
 * @brief This macro link the handler function to the SIGFPE. (Use this macro inside the main function)
 */
#define ATTACH_FPO_SIGNAL struct sigaction act = {};\
    act.sa_handler = InvalidFPOHandler;\
    sigaction(SIGFPE, &act, NULL);

#endif //DEBUG

#endif //WIN32

#endif //__FPO_EXCEPTIONS_H