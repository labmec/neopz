#ifndef CLOCK_TIMER_H

#include<string>
#include <config.h>
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
 
	int gettimeofday(struct timeval *tv, struct timezone *tz)
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
		tv->tv_sec = (long)(tmpres / 1000000UL);
		tv->tv_usec = (long)(tmpres % 1000000UL);
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
#else
#include<sys/time.h>
#endif

class ClockTimer
{
 public:
  ClockTimer() {};
  
  void start() {
    gettimeofday(&t1, NULL);
  }

  void stop() {
    timeval t2;
    gettimeofday(&t2, NULL);
    elapsed = (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
    elapsed += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms
  }

  double getUnits() {return (double) elapsed;}
    
  std::string getTime() {
    std::string str = double_to_string(elapsed); 
    str += " ms"; 
    return str;
  }

  private:

  timeval t1;
  double elapsed;
    
  std::string double_to_string(double dbl) {
    std::ostringstream strs;
    strs << dbl;
    return strs.str();
  }

};

#endif
