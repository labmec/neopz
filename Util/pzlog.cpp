/** 
 * @file 
 * @brief Contains the implementation of the InitializePZLOG() function. 
 */

#include "pzlog.h"
#include <sys/stat.h>
#include <iostream>

#ifdef LOG4CXX
std::mutex glogmutex;
#endif

void InitializePZLOG()
{
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
	int res = mkdir ("LOG", S_IRWXU | S_IXGRP | S_IRGRP | S_IXOTH | S_IROTH);
	// Wether the error happen again, the problem can to be permission, then a message is printed
	if(res) {
		struct stat result;
		res = stat ("LOG",&result);
		if(res) std::cout << "Error in mkdir : " << res << " permission denied." << std::endl;
	}
#endif
	
	std::cout << "Logfile " << configfile << std::endl;
	InitializePZLOG(configfile);
}

