/**
 * @file
 * @brief Defines PZError
 */
#ifndef PZERRORH
#define PZERRORH

#include <iostream>
#include <stdlib.h>

/**
 * @ingroup common
 * @brief Defines the output device to error messages and the DebugStop() function.
 */
#define PZError std::cout
/**
 * @ingroup common
 * @brief Returns a message to user put a breakpoint in
 */
#define DebugStop() pzinternal::DebugStopImpl(__FILE__, __LINE__)

namespace pzinternal {
void DebugStopImpl(const char *fileName, const std::size_t lineN);
}

#endif
