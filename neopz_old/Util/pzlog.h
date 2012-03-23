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

extern pthread_mutex_t glogmutex;

/*    pthread_mutex_lock(&fCommunicate);
    int ret = msg.ReceiveBlocking();
    pthread_mutex_unlock(&fCommunicate);
*/

#define LOGPZ_DEBUG(A,B) {pthread_mutex_lock(&glogmutex); \
                          LOG4CXX_DEBUG(A,B); \
                          pthread_mutex_unlock(&glogmutex); }
#define LOGPZ_INFO(A,B) {pthread_mutex_lock(&glogmutex); \
                          LOG4CXX_INFO(A,B) \
                          pthread_mutex_unlock(&glogmutex); }
#define LOGPZ_WARN(A,B) {pthread_mutex_lock(&glogmutex); \
                        LOG4CXX_WARN(A,B) \
                          pthread_mutex_unlock(&glogmutex); }
#define LOGPZ_ERROR(A,B) {pthread_mutex_lock(&glogmutex); \
                        LOG4CXX_ERROR(A,B) \
                          pthread_mutex_unlock(&glogmutex); }
#define LOGPZ_FATAL(A,B) {pthread_mutex_lock(&glogmutex); \
                        LOG4CXX_FATAL(A,B) \
                          pthread_mutex_unlock(&glogmutex); }

#else

#define LOGPZ_DEBUG(A,B) {}
#define LOGPZ_INFO(A,B) {}
#define LOGPZ_WARN(A,B) {}
#define LOGPZ_ERROR(A,B) {std::cout << B << std::endl;}
#define LOGPZ_FATAL(A,B) {std::cout << B << std::endl;}


#endif

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


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

void InitializePZLOG();

#endif
