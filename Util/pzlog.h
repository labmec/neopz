/**
 * @file pzlog.h
 * @brief Contains definitions of TPZLogger class, InitializePZLOG()  and auxiliary macros.

  TPZLogger is a wrapper for the logging library used in NeoPZ.
  The log config file is available at PZ_INSTALL_INCLUDE__DIR/Util/log4cxx.cfg, and the log output is at a folder LOG at **your executable** directory, not the directory of NeoPZ library.
  The auxiliary macros are LOGPZ_DEBUG, LOGPZ_INFO, LOGPZ_WARN, LOGPZ_ERROR and LOGPZ_FATAL.
  The functions TPZLogger::isXXXEnabled, with XXX=Info,Debug,Warn,Error and False should ALWAYS be called before writing the log messages, as most of the execution time is spent writing the messages to the stream.
 */


/**
 * \addtogroup util
 * \{ */

#ifndef PZLOGH
#define PZLOGH

#include <pz_config.h>


#ifdef PZ_LOG
#include <string>
class TPZLogger;
void LogPzDebugImpl(TPZLogger lg, std::string msg,
                    [[maybe_unused]] const char *fileName,[[maybe_unused]] const std::size_t lineN);
void LogPzInfoImpl(TPZLogger lg, std::string msg,
                   [[maybe_unused]] const char *fileName, [[maybe_unused]] const std::size_t lineN);
void LogPzWarnImpl(TPZLogger lg, std::string msg,
                   [[maybe_unused]] const char *fileName, [[maybe_unused]] const std::size_t lineN);
void LogPzErrorImpl(TPZLogger lg, std::string msg,
                    [[maybe_unused]] const char *fileName, [[maybe_unused]] const std::size_t lineN);
void LogPzFatalImpl(TPZLogger lg, std::string msg, 
                    [[maybe_unused]] const char *fileName, [[maybe_unused]] const std::size_t lineN);


/*The TPZLogger is a wrapper for the Log4cxx library. Once the
 log is setup, the expected usage is (in .cpp file):
 static TPZLogger logger("log_name");
 //do your stuff
 if(logger.isDebugEnabled()){
   //compose message
   LOGPZ_DEBUG(logger,message)
 }
 }*/

class TPZLogger{
public:
  TPZLogger() = delete;
  //get logger name as input param
  TPZLogger(const std::string &&);
  bool isDebugEnabled() const {return fIsDebugEnabled;}
  bool isWarnEnabled() const {return fIsWarnEnabled;}
  bool isInfoEnabled() const {return fIsInfoEnabled;}
  bool isErrorEnabled() const {return fIsErrorEnabled;}
  bool isFatalEnabled() const {return fIsFatalEnabled;}
private:
  const std::string fLogName;
  bool fIsDebugEnabled;
  bool fIsWarnEnabled;
  bool fIsInfoEnabled;
  bool fIsErrorEnabled;
  bool fIsFatalEnabled;

  friend void LogPzDebugImpl( TPZLogger lg,  std::string msg,
                     const char *fileName,  const std::size_t lineN);
  friend void LogPzInfoImpl( TPZLogger lg,  std::string msg,
                    const char *fileName,  const std::size_t lineN);
  friend void LogPzWarnImpl( TPZLogger lg,  std::string msg,
                    const char *fileName,  const std::size_t lineN);
  friend void LogPzErrorImpl( TPZLogger lg,  std::string msg,
                     const char *fileName,  const std::size_t lineN);
  friend void LogPzFatalImpl( TPZLogger lg,  std::string msg,
                     const char *fileName,  const std::size_t lineN);
};

/// Define log for debug
#define LOGPZ_DEBUG(A,B) LogPzDebugImpl(A,B,__FILE__,__LINE__);

/// Define log for info
#define LOGPZ_INFO(A,B) LogPzInfoImpl(A,B,__FILE__,__LINE__);

/// Define log for warnings
#define LOGPZ_WARN(A,B) LogPzWarnImpl(A,B,__FILE__,__LINE__);

/// Define log for errors
#define LOGPZ_ERROR(A,B) LogPzErrorImpl(A,B,__FILE__,__LINE__);
/// Define log for fatal errors
#define LOGPZ_FATAL(A,B) LogPzFatalImpl(A,B,__FILE__,__LINE__);

#else
#include <iostream>
//dummy class. log is not enabled.
class TPZLogger{
public:
  TPZLogger() = delete;
  
  inline TPZLogger(const std::string &&) {};
  inline bool isDebugEnabled() const {return false;}
  inline bool isWarnEnabled() const {return false;}
  inline bool isInfoEnabled() const {return false;}
  inline bool isErrorEnabled() const {return false;}
  inline bool isFatalEnabled() const {return false;}
};
// Just to allow the macros being called regardless
// of whether the log is enabled or not
#define LOGPZ_DEBUG(A,B)

#define LOGPZ_INFO(A,B)

#define LOGPZ_WARN(A,B)

#define LOGPZ_ERROR(A,B) \
std::cout<<B<<std::endl;\
DebugStop();

#define LOGPZ_FATAL(A,B)\
std::cout<<B<<std::endl;\
DebugStop();

#endif

/**
 * @brief Initializes log file for log4cxx with commom name log4cxx.cfg
 * @ingroup util
 */
void InitializePZLOG();

#endif
