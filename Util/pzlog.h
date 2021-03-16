/**
 * @file
 * @brief Contains definitions to LOGPZ_DEBUG, LOGPZ_INFO, LOGPZ_WARN, LOGPZ_ERROR and LOGPZ_FATAL,
 * and the implementation of the inline InitializePZLOG(string) function using log4cxx library or not.\n
 * It must to be called out of "#ifdef LOG4CXX" scope.
 */

#ifndef PZLOGH
#define PZLOGH

#include <string>
#include <sstream>
#include <iostream>           


#include <pz_config.h>


#ifdef LOG4CXX

#include <mutex>
#include <iostream>
#include <log4cxx/logger.h>
#include <log4cxx/basicconfigurator.h>
#include <log4cxx/propertyconfigurator.h>
using namespace log4cxx;
using namespace log4cxx::helpers;


/**
 * \addtogroup util
 * \{ */

/// External variable to mutex which controls write log
extern std::mutex glogmutex;

///EBORIN: PERF FIX: These macros  lock and unlock mutexes even if the
///        LOG4CXX macro does not have to log anything. Assuming the log
///        level does not change during execution, we could check for log level
///        before locking.
///PHIL : PERFORMANCE NOTE : most of the time is spent composing the string for logging
///        Check on the loglevel before composing the log string

/// Define log for debug
#define LOGPZ_DEBUG(A,B) {                      \
    if(A->isDebugEnabled()){                    \
      std::scoped_lock lock(glogmutex);         \
      LOG4CXX_DEBUG(A,B);                       \
    }                                           \
    else \
    {\
        std::cout << "Coloque IsDebugEnabled em " << __FILE__ << ":" << __LINE__ << std::endl;\
    }}

/// Define log for info
#define LOGPZ_INFO(A,B) {                                       \
    if(A->isInfoEnabled()) {                                    \
      std::scoped_lock lock(glogmutex);                         \
      LOG4CXX_INFO(A,B);}}

/// Define log for warnings
#define LOGPZ_WARN(A,B) {                                               \
    if(A->isWarnEnabled()) {                                            \
      std::scoped_lock lock(glogmutex);                                 \
      LOG4CXX_WARN(A,B);}}
/// Define log for errors

#define LOGPZ_ERROR(A,B) {                                      \
    if(A->isErrorEnabled()) {                                   \
      std::scoped_lock lock(glogmutex);                         \
      LOG4CXX_ERROR(A,B); DebugStop();}}

/// Define log for fatal errors
#define LOGPZ_FATAL(A,B) {                                  \
    if(A->isFatalEnabled()) {                               \
      std::scoped_lock lock(glogmutex);                     \
      LOG4CXX_FATAL(A,B);}}

#else

/// Define log for debug info
#define LOGPZ_DEBUG(A,B) {}
/// Define log for informations
#define LOGPZ_INFO(A,B) {}
/// Define log for warnings
#define LOGPZ_WARN(A,B) {}
/// Define log for errors (cout)
#define LOGPZ_ERROR(A,B) {std::cout << B << std::endl;}
/// Define log for fatal errors (cout)
#define LOGPZ_FATAL(A,B) {std::cout << B << std::endl;}


#endif


/**
 * @brief Initializes log file for log4cxx with commom name log4cxx.cfg
 * @ingroup util
 */
void InitializePZLOG();

#endif
