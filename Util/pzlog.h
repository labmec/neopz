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

inline void StopError()
{
	std::cout << "Ponto de parada\n";
}


#ifdef LOG4CXX

#include "pz_pthread.h"

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
extern pthread_mutex_t glogmutex;

///EBORIN: PERF FIX: These macros  lock and unlock mutexes even if the
///        LOG4CXX macro does not have to log anything. Assuming the log
///        level does not change during execution, we could check for log level
///        before locking.
///PHIL : PERFORMANCE NOTE : most of the time is spent composing the string for logging
///        Check on the loglevel before composing the log string

/// Define log for debug
#define LOGPZ_DEBUG(A,B) { \
      if(A->isDebugEnabled()) {		\
      PZ_PTHREAD_MUTEX_LOCK(&glogmutex,"LOGPZ_DEBUG");	   \
      LOG4CXX_DEBUG(A,B);				   \
      PZ_PTHREAD_MUTEX_UNLOCK(&glogmutex,"LOGPZ_DEBUG"); } \
        else std::cout << "Coloque IsDebugEnabled em " << __FILE__ << ":" << __LINE__ << std::endl;}

/// Define log for info
#define LOGPZ_INFO(A,B) {if(A->isInfoEnabled()) {		\
      PZ_PTHREAD_MUTEX_LOCK(&glogmutex, "LOGPZ_INFO");		\
      LOG4CXX_INFO(A,B);					\
      PZ_PTHREAD_MUTEX_UNLOCK(&glogmutex, "LOGPZ_INFO"); } }

/// Define log for warnings
#define LOGPZ_WARN(A,B) {if(A->isWarnEnabled()) {			\
      PZ_PTHREAD_MUTEX_LOCK(&glogmutex,"LOGPZ_WARN");			\
      LOG4CXX_WARN(A,B);						\
      PZ_PTHREAD_MUTEX_UNLOCK(&glogmutex,"LOGPZ_WARN"); } }
/// Define log for errors

#define LOGPZ_ERROR(A,B) {if(A->isErrorEnabled()) {		\
      PZ_PTHREAD_MUTEX_LOCK(&glogmutex,"LOGPZ_ERROR");		\
      LOG4CXX_ERROR(A,B); StopError();				\
      PZ_PTHREAD_MUTEX_UNLOCK(&glogmutex,"LOGPZ_ERROR"); } }

/// Define log for fatal errors
#define LOGPZ_FATAL(A,B) {if(A->isFatalEnabled()) {         \
      PZ_PTHREAD_MUTEX_LOCK(&glogmutex,"LOGPZ_FATAL");	    \
      LOG4CXX_FATAL(A,B);				    \
      PZ_PTHREAD_MUTEX_LOCK(&glogmutex,"LOGPZ_FATAL"); } }

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

/** \} */

/**
 * @ingroup util
 * @brief Initialize a log file adequated to use log4cxx lib
 */
inline void InitializePZLOG(const std::string &configfile)
{
#ifdef LOG4CXX
	log4cxx::PropertyConfigurator::configure(configfile);
	{
		log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("pz.mesh.tpzgeoelrefpattern"));
		logger->setAdditivity(false);
		//    logger->setLevel(log4cxx::Level::getDebug());
	}
	{
		log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("pz.mesh.refpattern"));
		logger->setAdditivity(false);
		//  logger->setLevel(log4cxx::Level::getDebug());
	}
#endif
}

/**
 * @brief Initializes log file for log4cxx with commom name log4cxx.cfg
 * @ingroup util
 */
void InitializePZLOG();

#endif
