/**
 * @file
 * @brief Contains definitions to LOGPZ_DEBUG, LOGPZ_INFO, LOGPZ_WARN, LOGPZ_ERROR and LOGPZ_FATAL, \n
 * and the implementation of the inline InitializePZLOG(string ) function.
 */
//
// C++ Interface: pzlog
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2005
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef PZLOGH
#define PZLOGH

#include <string>
#include <sstream>
#ifdef LOG4CXX

#include <pthread.h>

#ifdef DEBUG
  #define DEBUG2
#endif

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

/// Define log for debug
#define LOGPZ_DEBUG(A,B) {pthread_mutex_lock(&glogmutex); \
                          LOG4CXX_DEBUG(A,B); \
                          pthread_mutex_unlock(&glogmutex); }
/// Define log for info
#define LOGPZ_INFO(A,B) {pthread_mutex_lock(&glogmutex); \
                          LOG4CXX_INFO(A,B) \
                          pthread_mutex_unlock(&glogmutex); }
/// Define log for warnings
#define LOGPZ_WARN(A,B) {pthread_mutex_lock(&glogmutex); \
                        LOG4CXX_WARN(A,B) \
                          pthread_mutex_unlock(&glogmutex); }
/// Define log for errors
#define LOGPZ_ERROR(A,B) {pthread_mutex_lock(&glogmutex); \
                        LOG4CXX_ERROR(A,B) \
                          pthread_mutex_unlock(&glogmutex); }
/// Define log for fatal errors
#define LOGPZ_FATAL(A,B) {pthread_mutex_lock(&glogmutex); \
                        LOG4CXX_FATAL(A,B) \
                          pthread_mutex_unlock(&glogmutex); }

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

#ifdef HAVE_CONFIG_H
#include <config.h>
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
