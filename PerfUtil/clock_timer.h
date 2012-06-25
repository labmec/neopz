#ifndef CLOCK_TIMER_H

#include<pz_gettime.h>

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
