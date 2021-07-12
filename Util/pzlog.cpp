/** 
 * @file 
 * @brief Contains the implementation of the InitializePZLOG() function. 
 */

#include "pzlog.h"
#include "pzerror.h"
#include <sys/stat.h>
#include <iostream>



#ifdef PZ_LOG
#include <mutex>
static std::mutex glogmutex;

#include <log4cxx/logger.h>
#include <log4cxx/basicconfigurator.h>
#include <log4cxx/propertyconfigurator.h>




TPZLogger::TPZLogger(const std::string &&loggerName) : fLogName(loggerName){
}

void TPZLogger::InitializeLogLevels()
{
     auto logPtr = log4cxx::LoggerPtr(log4cxx::Logger::getLogger(fLogName));
     fIsDebugEnabled = logPtr->isDebugEnabled();
     fIsWarnEnabled = logPtr->isWarnEnabled();
     fIsInfoEnabled = logPtr->isInfoEnabled();
     fIsErrorEnabled = logPtr->isErrorEnabled();
     fIsFatalEnabled = logPtr->isFatalEnabled();
     fLogNotInitialized = false;
}


void pzinternal::LogPzDebugImpl(TPZLogger pzlg, std::string msg,
                                const char *funcName,const char *fileName,
                                const std::size_t lineN){
  log4cxx::LoggerPtr lg = log4cxx::LoggerPtr(log4cxx::Logger::getLogger(pzlg.fLogName));
  log4cxx::helpers::MessageBuffer oss_;
  lg->forcedLog(log4cxx::Level::getDebug(), oss_.str(oss_ << msg),
                log4cxx::spi::LocationInfo(fileName, funcName,lineN));
}

void pzinternal::LogPzInfoImpl(TPZLogger pzlg, std::string msg,
                               const char *funcName,const char *fileName,
                               const std::size_t lineN){
  log4cxx::LoggerPtr lg = log4cxx::LoggerPtr(log4cxx::Logger::getLogger(pzlg.fLogName));
  log4cxx::helpers::MessageBuffer oss_;
  lg->forcedLog(log4cxx::Level::getInfo(), oss_.str(oss_ << msg),
                log4cxx::spi::LocationInfo(fileName,funcName,lineN));
}
void pzinternal::LogPzWarnImpl(TPZLogger pzlg, std::string msg,
                               const char *funcName, const char *fileName,
                               const std::size_t lineN){
  log4cxx::LoggerPtr lg = log4cxx::LoggerPtr(log4cxx::Logger::getLogger(pzlg.fLogName));
  log4cxx::helpers::MessageBuffer oss_;
  lg->forcedLog(log4cxx::Level::getWarn(), oss_.str(oss_ << msg),
                log4cxx::spi::LocationInfo(fileName, funcName,lineN));
}

void pzinternal::LogPzErrorImpl(TPZLogger pzlg, std::string msg,
                                const char *funcName, const char *fileName,
                                const std::size_t lineN){
  log4cxx::LoggerPtr lg = log4cxx::LoggerPtr(log4cxx::Logger::getLogger(pzlg.fLogName));
  log4cxx::helpers::MessageBuffer oss_;
  lg->forcedLog(log4cxx::Level::getError(), oss_.str(oss_ << msg),
                log4cxx::spi::LocationInfo(fileName,funcName,lineN));
}

void pzinternal::LogPzFatalImpl(TPZLogger pzlg, std::string msg,
                                const char *funcName, const char *fileName,
                                const std::size_t lineN) {
  log4cxx::LoggerPtr lg = log4cxx::LoggerPtr(log4cxx::Logger::getLogger(pzlg.fLogName));
  log4cxx::helpers::MessageBuffer oss_; 
  lg->forcedLog(::log4cxx::Level::getFatal(), oss_.str(oss_ << msg),
                log4cxx::spi::LocationInfo(fileName,funcName,lineN));
}

/**
 * @ingroup util
 * @brief Initialize a log file adequated to use log4cxx lib
 */
void TPZLogger::InitializePZLOG(const std::string &configfile)
{
  std::cout << "Using the following NeoPZ log config file:" << std::endl;
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
  std::cout  << configfile << std::endl;

  log4cxx::PropertyConfigurator::configure(configfile);
}

#endif
