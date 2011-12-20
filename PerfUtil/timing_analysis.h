#ifndef TIMING_ANALYSIS_H

#include<string>
#include<set>
#include<iostream>

#define TIME_ANALYSIS_ADD(ta, timer, ...)                        \
  { std::ostringstream ostr;                                     \
    ostr << __VA_ARGS__;					 \
    ta.add(timer.getUnits(), ostr.str()); }

#define TIME_SEC_BEG_LOG(logger,timer, ...)		   \
  LOG4CXX_INFO(logger, __VA_ARGS__ << " started")	   \
  timer.start()

#define TIME_SEC_END_LOG(logger,ta,timer, ...)				\
  timer.stop();                                                         \
  TIME_ANALYSIS_ADD(ta,timer, __VA_ARGS__)                              \
  LOG4CXX_INFO(logger, __VA_ARGS__ << " ended in " << timer.getTime())

#define TIME_SEC_BEG(timer, ...)                   \
  timer.start()

#define TIME_SEC_END(ta,timer, ...)                                     \
  timer.stop();                                                         \
  TIME_ANALYSIS_ADD(ta,timer, __VA_ARGS__)

typedef std::pair<double,std::string> section_timing_t;
typedef std::set<section_timing_t> timing_set_t;

class TimingAnalysis
{
 public:
  TimingAnalysis() {};
  
  void add(double units, const std::string& section_name) {
    section_timing_t st(units,section_name);
    ts.insert(st);
  }

  void share_report(std::ostream& o, double total) {
    double share;
    double sum=0;
    o << "Timing analysis Report" << std::endl;
    timing_set_t::iterator it;
    for(it = ts.begin(); it != ts.end(); it++) {
      sum += it->first;
      share = 100 * (it->first / total);
      o << std::setw(6) << std::setprecision(2)
	<< std::fixed << share << " % : ";
      o << it->second << " : " << it->first << std::endl; 
    }
    share = 100 * (sum / total);
    o << std::setw(6) << std::setprecision(2)
      << std::fixed << share << " % : Total";
    o << " : " << sum << std::endl; 
  }

  void share_report(std::ostream& o) {
    double total = 0;
    timing_set_t::iterator it;
    for(it = ts.begin(); it != ts.end(); it++) {
      total += it->first;
    }
    share_report(o, total);
  }

  private:

  timing_set_t ts;
 
};


#endif
