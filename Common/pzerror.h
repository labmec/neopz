#ifndef PZERRORH
#define PZERRORH

#include <iostream>
#include <stdlib.h>

/** 
 * @ingroup common
 * @brief Defines the output device to error messages. (Jorge)
 */
#define PZError std::cout

/**
 * @ingroup common
 * @brief Returns a message to user put a breakpoint in 
 */
void DebugStop();

#endif
