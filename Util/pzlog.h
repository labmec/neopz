/**
 * @file pzlog.h
 * @brief Contains definitions of TPZLogger class, InitializePZLOG()  and auxiliary macros.

  The auxiliary macros are LOGPZ_DEBUG, LOGPZ_INFO, LOGPZ_WARN, LOGPZ_ERROR and LOGPZ_FATAL. The functions
  TPZLogger::isXXXEnabled, with XXX=Info,Debug,Warn,Error and False should ALWAYS be called before writing the log messages, as most of the time is spent composing the messages.
 */


/**
 * \addtogroup util
 * \{ */

#ifndef PZLOGH
#define PZLOGH

#include <string>
#include <sstream>
#include <iostream>           


#include <pz_config.h>
#include <mutex>
#include <iostream>




/// External variable to mutex which controls write log
#ifdef PZ_LOG
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



/*@TODO:Check if we cannot call getlogger in the functions called by
 the macros. If there are no considerable performance issues,
 then we could eliminate this forward declaration and the attribute
 in the TPZLogger class, storing just the name of the logger instead.*/
namespace log4cxx{
  namespace helpers{
    template<class T>
    class ObjectPtrT;
  }
  class Logger;
  typedef log4cxx::helpers::ObjectPtrT<Logger> LoggerPtr;
}



class TPZLogger{
public:
  TPZLogger() = delete;
  //get logger name as input param
  TPZLogger(const std::string &&);
  ~TPZLogger();
  bool isDebugEnabled() const;
  bool isWarnEnabled() const;
  bool isInfoEnabled() const;
  bool isErrorEnabled() const;
  bool isFatalEnabled() const;
private:
  log4cxx::LoggerPtr * fLogPtr;

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
extern std::mutex glogmutex;
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
