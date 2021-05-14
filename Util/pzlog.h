/**
 * @file pzlog.h
 * @brief Contains definitions of TPZLogger class, InitializePZLOG()  and auxiliary macros.
 *
 * TPZLogger is a wrapper for the logging library used in NeoPZ.
 * The default log output is at a folder LOG at **your executable** directory, not the directory of NeoPZ library.
 * Once the NeoPZ library is installed, the default log config file to be used is the one in PZ_INSTALL_DIR/include/Util/log4cxx.cfg.
 * The logger can be initialized with any custom config through the 
 * function TPZLogger:InitializePZLOG(const std::string&) receiving
 * the name (and path, obviously) of the config file.
 * The auxiliary macros are LOGPZ_DEBUG, LOGPZ_INFO, LOGPZ_WARN, LOGPZ_ERROR and LOGPZ_FATAL.
 * The functions TPZLogger::isXXXEnabled, with XXX=Info,Debug,Warn,Error and False should ALWAYS be called
 * before writing the log messages, as most of the execution time is spent writing the messages to the stream.
 * If you are extremely worried about performance, you can use 
 * 
 * #ifdef PZ_LOG
 * //your log code here
 * #endif
 * 
 * on the log calls. However, it is really not necessary.
 * Most compilers are probably able to optimise the if(logger.isXXXenabled()) condition when PZ_LOG is not set.
 */




#ifndef PZLOGH
#define PZLOGH

#include <pz_config.h>

/**
 * \addtogroup util
 * \{ */

#ifdef PZ_LOG
#include <string>

class TPZLogger;

//the following are functions for internal usage

namespace pzinternal{
void LogPzDebugImpl(TPZLogger lg, std::string msg,
                    [[maybe_unused]] const char *fileName,
                    [[maybe_unused]] const std::size_t lineN);
void LogPzInfoImpl(TPZLogger lg, std::string msg,
                   [[maybe_unused]] const char *fileName,
                   [[maybe_unused]] const std::size_t lineN);
void LogPzWarnImpl(TPZLogger lg, std::string msg,
                   [[maybe_unused]] const char *fileName,
                   [[maybe_unused]] const std::size_t lineN);
void LogPzErrorImpl(TPZLogger lg, std::string msg,
                    [[maybe_unused]] const char *fileName,
                    [[maybe_unused]] const std::size_t lineN);
void LogPzFatalImpl(TPZLogger lg, std::string msg,
                    [[maybe_unused]] const char *fileName,
                    [[maybe_unused]] const std::size_t lineN);
}

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
  //the following functions are getters for logger lvls
  bool isDebugEnabled() const {return fIsDebugEnabled;}
  bool isWarnEnabled() const {return fIsWarnEnabled;}
  bool isInfoEnabled() const {return fIsInfoEnabled;}
  bool isErrorEnabled() const {return fIsErrorEnabled;}
  bool isFatalEnabled() const {return fIsFatalEnabled;}
  //initializes loggger using custom config in configfile.
  static void InitializePZLOG(const std::string &configfile);
  //initializes logger using default NeoPZ config
  static void InitializePZLOG(){
    std::string path = PZ_LOG4CXX_CONFIG_FILE;
    TPZLogger::InitializePZLOG(path);
  }
private:
  const std::string fLogName;
  bool fIsDebugEnabled;
  bool fIsWarnEnabled;
  bool fIsInfoEnabled;
  bool fIsErrorEnabled;
  bool fIsFatalEnabled;

  friend void pzinternal::LogPzDebugImpl(TPZLogger lg, std::string msg,
                                         const char *fileName,
                                         const std::size_t lineN);
  friend void pzinternal::LogPzInfoImpl(TPZLogger lg, std::string msg,
                                        const char *fileName,
                                        const std::size_t lineN);
  friend void pzinternal::LogPzWarnImpl(TPZLogger lg, std::string msg,
                                        const char *fileName,
                                        const std::size_t lineN);
  friend void pzinternal::LogPzErrorImpl(TPZLogger lg, std::string msg,
                                         const char *fileName,
                                         const std::size_t lineN);
  friend void pzinternal::LogPzFatalImpl(TPZLogger lg, std::string msg,
                                         const char *fileName,
                                         const std::size_t lineN);
};

//the following macros are for internal usage

/// Define log for debug
#define LOGPZ_DEBUG(A,B) pzinternal::LogPzDebugImpl(A,B,__FILE__,__LINE__);

/// Define log for info
#define LOGPZ_INFO(A,B) pzinternal::LogPzInfoImpl(A,B,__FILE__,__LINE__);

/// Define log for warnings
#define LOGPZ_WARN(A,B) pzinternal::LogPzWarnImpl(A,B,__FILE__,__LINE__);

/// Define log for errors
#define LOGPZ_ERROR(A,B) pzinternal::LogPzErrorImpl(A,B,__FILE__,__LINE__);
/// Define log for fatal errors
#define LOGPZ_FATAL(A,B) pzinternal::LogPzFatalImpl(A,B,__FILE__,__LINE__);


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
  //initializes loggger using custom config in configfile.
  static void InitializePZLOG(const std::string &configfile){}
  //initializes logger using default NeoPZ config
  static void InitializePZLOG(){}
};
// Just to allow the macros being called regardless
// of whether the log is enabled or not
#define LOGPZ_DEBUG(A,B)

#define LOGPZ_INFO(A,B)

#define LOGPZ_WARN(A,B)

#define LOGPZ_ERROR(A,B)   \
{                          \
  std::cout<<B<<std::endl; \
  DebugStop();             \
}

#define LOGPZ_FATAL(A,B)   \
{                          \
  std::cout<<B<<std::endl; \
  DebugStop();             \
}

#endif
/**
 * \} */
#endif
