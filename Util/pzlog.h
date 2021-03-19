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


#include <mutex>
#include <iostream>




/// External variable to mutex which controls write log
#ifdef PZ_LOG
class PZLogger;
void LogPzDebugImpl([[maybe_unused]] PZLogger lg, [[maybe_unused]] std::string msg,
                    [[maybe_unused]] const char *fileName, [[maybe_unused]] const std::size_t lineN);
void LogPzInfoImpl([[maybe_unused]] PZLogger lg, [[maybe_unused]] std::string msg,
                   [[maybe_unused]] const char *fileName, [[maybe_unused]] const std::size_t lineN);
void LogPzWarnImpl([[maybe_unused]] PZLogger lg, [[maybe_unused]] std::string msg,
                   [[maybe_unused]] const char *fileName, [[maybe_unused]] const std::size_t lineN);
void LogPzErrorImpl([[maybe_unused]] PZLogger lg, [[maybe_unused]] std::string msg,
                    [[maybe_unused]] const char *fileName, [[maybe_unused]] const std::size_t lineN);
void LogPzFatalImpl([[maybe_unused]] PZLogger lg, [[maybe_unused]] std::string msg,
                    [[maybe_unused]] const char *fileName, [[maybe_unused]] const std::size_t lineN);

namespace log4cxx{
  namespace helpers{
    template<class T>
    class ObjectPtrT;
  }
  class Logger;
  typedef log4cxx::helpers::ObjectPtrT<Logger> LoggerPtr;
}



class PZLogger{
public:
  PZLogger() = delete;
  
  PZLogger(const std::string &&);
  ~PZLogger();
  bool isDebugEnabled() const;
  bool isWarnEnabled() const;
  bool isInfoEnabled() const;
private:
  log4cxx::LoggerPtr * fLogPtr;

  friend void LogPzDebugImpl( PZLogger lg,  std::string msg,
                     const char *fileName,  const std::size_t lineN);
  friend void LogPzInfoImpl( PZLogger lg,  std::string msg,
                    const char *fileName,  const std::size_t lineN);
  friend void LogPzWarnImpl( PZLogger lg,  std::string msg,
                    const char *fileName,  const std::size_t lineN);
  friend void LogPzErrorImpl( PZLogger lg,  std::string msg,
                     const char *fileName,  const std::size_t lineN);
  friend void LogPzFatalImpl( PZLogger lg,  std::string msg,
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
class PZLogger{
public:
  PZLogger() = delete;
  
  inline PZLogger(const std::string &&) {};
  bool isDebugEnabled() const {return false;}
  bool isWarnEnabled() const {return false;}
  bool isInfoEnabled() const {return false;}
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
