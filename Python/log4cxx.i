// SWIG interface
%module log4cxx
%{
// #include "/usr/local/include/log4cxx/helpers/object.h"
// #include "/usr/local/include/log4cxx/log4cxx.h"
// #include "/usr/local/include/log4cxx/logger.h"
// #include "/usr/local/include/log4cxx/helpers/appenderattachableimpl.h"
// #include "/usr/local/include/log4cxx/level.h"
// #include "/usr/local/include/log4cxx/helpers/pool.h"
// #include "/usr/local/include/log4cxx/helpers/mutex.h"
// #include "/usr/local/include/log4cxx/spi/location/locationinfo.h"
// #include "/usr/local/include/log4cxx/helpers/resourcebundle.h"
// #include "/usr/local/include/log4cxx/helpers/messagebuffer.h"

#include <log4cxx/helpers/object.h>
#include <log4cxx/log4cxx.h>
#include <log4cxx/logger.h>
#include <log4cxx/helpers/appenderattachableimpl.h>
#include <log4cxx/level.h>
#include <log4cxx/helpers/pool.h>
#include <log4cxx/helpers/mutex.h>
#include <log4cxx/spi/location/locationinfo.h>
#include <log4cxx/helpers/resourcebundle.h>
#include <log4cxx/helpers/messagebuffer.h>

#include <log4cxx/logger.h>
#include <log4cxx/basicconfigurator.h>
#include <log4cxx/propertyconfigurator.h>
using namespace log4cxx;
using namespace log4cxx::helpers;

LoggerPtr logger(log4cxx::Logger::getLogger("pz.mesh.tpzgeoelrefpattern"));
logger->setAdditivity(false);
%}

//%include "/usr/local/include/log4cxx/log4cxx.h"
//%include "/usr/local/include/log4cxx/logger.h"

