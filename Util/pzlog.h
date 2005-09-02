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

#ifdef LOG4CXX

#include "log4cxx/logger.h"
using namespace log4cxx;
using namespace log4cxx::helpers;

#define LOGPZ_DEBUG(A,B) LOG4CXX_DEBUG(A,B)
#define LOGPZ_INFO(A,B) LOG4CXX_INFO(A,B)
#define LOGPZ_WARN(A,B) LOG4CXX_WARN(A,B)
#define LOGPZ_ERROR(A,B) LOG4CXX_ERROR(A,B)
#define LOGPZ_FATAL(A,B) LOG4CXX_FATAL(A,B)



#else

#define LOGPZ_DEBUG(A,B) {}
#define LOGPZ_INFO(A,B) {}
#define LOGPZ_WARN(A,B) {}
#define LOGPZ_ERROR(A,B) {}
#define LOGPZ_FATAL(A,B) {}


#endif

