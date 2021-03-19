/** 
 * @file 
 * @brief Contains the implementation of the InitializePZLOG() function. 
 */

#include "pzlog.h"
#include "pzerror.h"
#include <sys/stat.h>
#include <iostream>


std::mutex glogmutex;

#ifdef PZ_LOG


#include <log4cxx/logger.h>
#include <log4cxx/basicconfigurator.h>
#include <log4cxx/propertyconfigurator.h>




TPZLogger::TPZLogger(const std::string &&loggerName) : fLogName(loggerName){  
   auto logPtr = log4cxx::LoggerPtr(log4cxx::Logger::getLogger(loggerName));
   fIsDebugEnabled = logPtr->isDebugEnabled();
   fIsWarnEnabled = logPtr->isWarnEnabled();
   fIsInfoEnabled = logPtr->isInfoEnabled();
   fIsErrorEnabled = logPtr->isErrorEnabled();
   fIsFatalEnabled = logPtr->isFatalEnabled();
}

void LogPzDebugImpl(TPZLogger pzlg, std::string msg, const char *fileName, const std::size_t lineN){
  log4cxx::LoggerPtr lg = log4cxx::LoggerPtr(log4cxx::Logger::getLogger(pzlg.fLogName));
  if (lg->isDebugEnabled()) {
    std::scoped_lock lock(glogmutex);
    LOG4CXX_DEBUG(lg,msg);
  } else {
    std::cout << "Please set isDebugEnabled at " << fileName << ":" << lineN
              << std::endl;
  }
}

void LogPzInfoImpl(TPZLogger pzlg, std::string msg, const char *fileName, const std::size_t lineN){
  log4cxx::LoggerPtr lg = log4cxx::LoggerPtr(log4cxx::Logger::getLogger(pzlg.fLogName));
  if (lg->isInfoEnabled()) {
    std::scoped_lock lock(glogmutex);
    LOG4CXX_INFO(lg, msg);
  } else {
    std::cout << "Please set isInfoEnabled at " << fileName << ":" << lineN
              << std::endl;
  }
}
void LogPzWarnImpl(TPZLogger pzlg, std::string msg, const char *fileName, const std::size_t lineN){
  log4cxx::LoggerPtr lg = log4cxx::LoggerPtr(log4cxx::Logger::getLogger(pzlg.fLogName));
  if (lg->isWarnEnabled()) {
    std::scoped_lock lock(glogmutex);
    LOG4CXX_WARN(lg, msg);
  }else {
    std::cout << "Please set isWarnEnabled at " << fileName << ":" << lineN
              << std::endl;
  }
}

void LogPzErrorImpl(TPZLogger pzlg, std::string msg, const char *fileName, const std::size_t lineN){
  log4cxx::LoggerPtr lg = log4cxx::LoggerPtr(log4cxx::Logger::getLogger(pzlg.fLogName));
  if (lg->isErrorEnabled()) {
    std::scoped_lock lock(glogmutex);
    LOG4CXX_ERROR(lg, msg);
    DebugStop();
  }else {
    std::cout << "Please set isErrorEnabled at " << fileName << ":" << lineN
              << std::endl;
  }
}

void LogPzFatalImpl(TPZLogger pzlg, std::string msg, const char *fileName, const std::size_t lineN){
  log4cxx::LoggerPtr lg = log4cxx::LoggerPtr(log4cxx::Logger::getLogger(pzlg.fLogName));
  if (lg->isFatalEnabled()) {
    std::scoped_lock lock(glogmutex);
    LOG4CXX_FATAL(lg, msg);
  }else {
    std::cout << "Please set isFatalEnabled at " << fileName << ":" << lineN
              << std::endl;
  }
}

#endif


/**
 * @ingroup util
 * @brief Initialize a log file adequated to use log4cxx lib
 */
inline void InitializePZLOG(const std::string &configfile)
{
#ifdef PZ_LOG
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

void InitializePZLOG()
{
  static bool firstTime{true};
  if(!firstTime){
    std::cout << "Are you calling InitializePZLOG()?"<< std::endl;
    std::cout << "This is not needed anymore"<<std::endl;
  }
  else{
    std::string path;
    std::string configfile;
#ifdef PZSOURCEDIR
    path = PZSOURCEDIR;
    path += "/Util/";
#else
    path = "";
#endif
    configfile = path;
    configfile += "log4cxx.cfg";

#ifndef WIN32
    int res = mkdir("LOG", S_IRWXU | S_IXGRP | S_IRGRP | S_IXOTH | S_IROTH);
    // Wether the error happen again, the problem can to be permission, then a
    // message is printed
    if (res) {
      struct stat result;
      res = stat("LOG", &result);
      if (res)
        std::cout << "Error in mkdir : " << res << " permission denied."
                  << std::endl;
    }
#endif

    std::cout << "Logfile " << configfile << std::endl;
    InitializePZLOG(configfile);
    firstTime = false;
  }
	
}


/**@orlandini: The following is an attempt to call
InitializePZLOG() at some point before main() is called.
If anybody comes up with a better solution, it would be great.
Taken from: https://stackoverflow.com/questions/19227664/whats-the-c-idiom-equivalent-to-the-java-static-block

Apparently, constructors of namespace-scope variables
are guaranteed to run before main. It could be
an alternative... See:
https://stackoverflow.com/questions/9439871/can-you-print-anything-in-c-before-entering-into-the-main-function*/


#define M_CON(A, B) M_CON_(A, B)
#define M_CON_(A, B) A##B

#define STATIC_BLOCK \
        [[maybe_unused]] static const auto M_CON(_static_block,__LINE__) = []()

STATIC_BLOCK {
    InitializePZLOG();
    return 0;
}();

#undef M_CON
#undef M_CON_
#undef STATIC_BLOCK
