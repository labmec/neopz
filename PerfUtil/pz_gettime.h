#ifndef PZ_GETTIME_H
#define PZ_GETTIME_H

#include <pz_config.h>

#ifdef VC
	//We need to implement gettimeofday function on windows environment.
	#include <time.h>
	#if defined(_MSC_VER) || defined(_MSC_EXTENSIONS)
	  #define DELTA_EPOCH_IN_MICROSECS  11644473600000000Ui64
	#else
	  #define DELTA_EPOCH_IN_MICROSECS  11644473600000000ULL
	#endif
 
	struct timezone 
	{
	  int  tz_minuteswest; /* minutes W of Greenwich */
	  int  tz_dsttime;     /* type of dst correction */
	};
 
	inline int gettimeofday(struct timeval *tv, struct timezone *tz)
	{
	  FILETIME ft;
	  unsigned __int64 tmpres = 0;
	  static int tzflag;
 
	  if (NULL != tv)
	  {
		GetSystemTimeAsFileTime(&ft);
 
		tmpres |= ft.dwHighDateTime;
		tmpres <<= 32;
		tmpres |= ft.dwLowDateTime;
 
		/*converting file time to unix epoch*/
		tmpres -= DELTA_EPOCH_IN_MICROSECS; 
		tmpres /= 10;  /*convert into microseconds*/
		tv->tv_sec = (int64_t)(tmpres / 1000000UL);
		tv->tv_usec = (int64_t)(tmpres % 1000000UL);
	  }
 
	  if (NULL != tz)
	  {
		if (!tzflag)
		{
		  _tzset();
		  tzflag++;
		}
		tz->tz_minuteswest = _timezone / 60;
		tz->tz_dsttime = _daylight;
	  }
 
	  return 0;
	}
#else // VC

#include<sys/time.h>

#endif // VC

#endif // PZ_GETTIME_H
