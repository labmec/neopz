/**
 * @file
 * @brief Contains macros and functions to support execution statistics recording.
 */
#ifndef STATS_RECORDER_H

#define HAS_GETRUSAGE

#include<sstream>  // stringstream
#include<string.h> // memset

/**
   Stats 
template<T>
class stat_table_col {

  const string getHeader() const;
  void setHeader(string header);
  
  // Return the number of rows
  unsigned nRows() const;
  
  // Add nrows rows to the column.
  unsigned addRows(unsigned nrows);

  const T& getCellValue(unsigned row) const;

  T& getCellValue(unsigned row);

  protected:

    string header;
    vector<T> data;
    
};
class stat_table {

   // Add a new row to the end and return the row id.
   unsigned addRow(); 

   // Set the value of a table cell. Return zero if ok, < 0 if error.
   // Error code: -1 : Invalid row
   //             -2 : Invalid column and createNewCol == false
   int setCell(unsigned row, string colId, string value, bool createNewCol = false);
   int setCell(unsigned row, string colId, int value, bool createNewCol = false);
   int setCell(unsigned row, string colId, bool value, bool createNewCol = false);
   int setCell(unsigned row, string colId, unsigned value, bool createNewCol = false);
   int setCell(unsigned row, string colId, double value, bool createNewCol = false);
   int setCell(unsigned row, string colId, float value, bool createNewCol = false);
  
   stat_table_col& getColumn(const string colId); 

   const stat_table_col& getColumn(const string colId) const; 
   
}

clarg::argString st_filename stfn("-stats", "Statistics table file name", "stats.csv");
stat_table st;

st.readFromFile(stfn.getValue());

  /* Contains statistics 
PerfAnalyzer pa_create_matrix;
PerfAnalyzer pa_assemble;

pa_create_matrix.start();
... // Profiled work
pa_create_matrix.stop();

... // Other stuff

pa_create_matrix.start();
... // Profiled work
pa_create_matrix.stop();

pa_create_matrix.
.n_laps;     //N_LAPS
.usr_time;   //USR_TIME
.sys_time;   //SYS_TIME
.real_time;  //REAL_TIME

- Append metrics to the statistics table.
- The idea is to keep one table for each segment of
  code. New runs are appended to the table.
pa_create_matrix.append_to(st);
   // unsigned new_row = st.addRow()
   // st.setCell(new_row, std::string("USR_TIME"), usr_time, true); 

st.writeToFile(stfn.getValue());

PerfAnalysisAggregate aggregated_pa;
aggregated_pa.add(pa_init);

PerfAnalysisTracker pat;

*/

using namespace std;

#include<csvtable.h>

/**
 * Base class to store statistics that should be diffed (eg. execution time)
 */
class RunStat
{
 public:
  /* Start recording the execution statistics. */
  virtual void start() = 0;
  /* Stop recording the execution statistics and accumulate the partial result. */
  virtual void stop() = 0;
  /**
   * Print the statistics  
   */
  virtual void print(ostream& os) const = 0;
  /**
   * Append metrics to the statistics table.  The idea is to keep one
   * table for each segment of code. New runs are appended to the
   * table.  
   * Returns: the new_row number if Ok 
   *          the setCell error code (< 0) if Error.
  */
  virtual int setCellValues(CSVStringTable& st, unsigned row) const = 0;

  /**
   * Reset the statistics. 
   */
  virtual void clearStats() = 0;
};

#ifdef HAS_GETRUSAGE
#include <sys/resource.h> // getrusage

class RUsageRunStat : public RunStat
{
 public:

  RUsageRunStat()
    {
      clearStats();
    }

  /* Start recording the execution statistics. */
  void start()
  {
    getrusage(RUSAGE_SELF, &lap_self);
    getrusage(RUSAGE_CHILDREN, &lap_children);
  }

  /* Stop recording the execution statistics. */
  void stop()
  {
    struct rusage self, children;
    getrusage(RUSAGE_SELF, &self);
    getrusage(RUSAGE_CHILDREN, &children);
#define SET_TOTAL(fld) total_self.fld += (self.fld - lap_self.fld); total_children.fld += (children.fld - lap_children.fld)
    //struct timeval ru_utime; /* user time used */
    //struct timeval ru_stime; /* system time used */
    SET_TOTAL(ru_maxrss);          /* integral max resident set size */
    SET_TOTAL(ru_ixrss);           /* integral shared text memory size */
    SET_TOTAL(ru_idrss);           /* integral unshared data size */
    SET_TOTAL(ru_isrss);           /* integral unshared stack size */
    SET_TOTAL(ru_minflt);          /* page reclaims */
    SET_TOTAL(ru_majflt);          /* page faults */
    SET_TOTAL(ru_nswap);           /* swaps */
    SET_TOTAL(ru_inblock);         /* block input operations */
    SET_TOTAL(ru_oublock);         /* block output operations */
    SET_TOTAL(ru_msgsnd);          /* messages sent */
    SET_TOTAL(ru_msgrcv);          /* messages received */
    SET_TOTAL(ru_nsignals);        /* signals received */
    SET_TOTAL(ru_nvcsw);           /* voluntary context switches */
    SET_TOTAL(ru_nivcsw);          /* involuntary context switches */
  }

  /**
   * Print the statistics.
   * TODO: create a table, update it, and print the table.
   */
  void print(ostream& os) const
  {
    stringstream header;
    stringstream values;
#define PRINT_FLD(hd,fld) header << ",SELF_" << hd << ",CHD_" << hd; values << "," << total_self.fld << "," << total_children.fld

    //struct timeval ru_utime; /* user time used */
    //struct timeval ru_stime; /* system time used */
    PRINT_FLD("RU_MAXRSS",ru_maxrss);          /* integral max resident set size */
    PRINT_FLD("RU_IXRSS",ru_ixrss);           /* integral shared text memory size */
    PRINT_FLD("RU_IDRSS",ru_idrss);           /* integral unshared data size */
    PRINT_FLD("RU_ISRSS",ru_isrss);           /* integral unshared stack size */
    PRINT_FLD("RU_MINFLT",ru_minflt);          /* page reclaims */
    PRINT_FLD("RU_MAJFLT",ru_majflt);          /* page faults */
    PRINT_FLD("RU_NSWAP",ru_nswap);           /* swaps */
    PRINT_FLD("RU_INBLOCK",ru_inblock);         /* block input operations */
    PRINT_FLD("RU_OUBLOCK",ru_oublock);         /* block output operations */
    PRINT_FLD("RU_MSGND",ru_msgsnd);          /* messages sent */
    PRINT_FLD("RU_MSGRCV",ru_msgrcv);          /* messages received */
    PRINT_FLD("RU_NSIGNAL",ru_nsignals);        /* signals received */
    PRINT_FLD("RU_NVCSW",ru_nvcsw);           /* voluntary context switches */
    PRINT_FLD("RU_NIVSW",ru_nivcsw);          /* involuntary context switches */
    os << "HEADERS" << header.str() << endl;
    os << "VALUES" << values.str() << endl;
  }

  /**
   * Append metrics to the statistics table.  The idea is to keep one
   * table for each segment of code. New runs are appended to the
   * table.  
   * Returns: the new_row number if Ok 
   *          the setCell error code (< 0) if Error.
  */
  int setCellValues(CSVStringTable& st, unsigned row) const 
  {
    if (st.nRows() <= row) return -1;

#define APPEND_LONG_FLD(hd,fld)						\
    { int ret; string header("SELF_"); header += hd;			\
      if ((ret = st.setCell(row, header, total_self.fld, true)))	\
	return ret;							\
      header = "CHD_"; header += hd;					\
      if ((ret = st.setCell(row, header, total_children.fld, true)))	\
	return ret;							\
    }

    //struct timeval ru_utime; /* user time used */
    //struct timeval ru_stime; /* system time used */
    APPEND_LONG_FLD("RU_MAXRSS",ru_maxrss);          /* integral max resident set size */
    APPEND_LONG_FLD("RU_IXRSS",ru_ixrss);           /* integral shared text memory size */
    APPEND_LONG_FLD("RU_IDRSS",ru_idrss);           /* integral unshared data size */
    APPEND_LONG_FLD("RU_ISRSS",ru_isrss);           /* integral unshared stack size */
    APPEND_LONG_FLD("RU_MINFLT",ru_minflt);          /* page reclaims */
    APPEND_LONG_FLD("RU_MAJFLT",ru_majflt);          /* page faults */
    APPEND_LONG_FLD("RU_NSWAP",ru_nswap);           /* swaps */
    APPEND_LONG_FLD("RU_INBLOCK",ru_inblock);         /* block input operations */
    APPEND_LONG_FLD("RU_OUBLOCK",ru_oublock);         /* block output operations */
    APPEND_LONG_FLD("RU_MSGND",ru_msgsnd);          /* messages sent */
    APPEND_LONG_FLD("RU_MSGRCV",ru_msgrcv);          /* messages received */
    APPEND_LONG_FLD("RU_NSIGNAL",ru_nsignals);        /* signals received */
    APPEND_LONG_FLD("RU_NVCSW",ru_nvcsw);           /* voluntary context switches */
    APPEND_LONG_FLD("RU_NIVSW",ru_nivcsw);          /* involuntary context switches */

    return 0; // Return ok
  }

  void clearStats()
  {
    memset(&lap_self, 0, sizeof(lap_self));
    memset(&lap_children, 0, sizeof(lap_children));
    memset(&total_self, 0, sizeof(total_self));
    memset(&total_children, 0, sizeof(total_children));
  }

 protected:

  struct rusage lap_self, lap_children;
  struct rusage total_self, total_children;

};

#endif

#include<pz_gettime.h>

class ElapsedTimeRunStat : public RunStat
{
 public:

  ElapsedTimeRunStat()
    {
      clearStats();
    }

  /* Start recording the execution statistics. */
  void start()
  {
    gettimeofday(&lap, NULL);
  }

  /* Stop recording the execution statistics. */
  void stop()
  {
    timeval curr;
    gettimeofday(&curr, NULL);
    total.tv_sec  += (curr.tv_sec  - lap.tv_sec );
    total.tv_usec += (curr.tv_usec - lap.tv_usec);
    //elapsed = (t.tv_sec - lap.tv_sec) * 1000.0;      // sec to ms
    //elapsed += (t.tv_usec - lap.tv_usec) / 1000.0;   // us to ms
  }

  /**
   * Print the statistics  
   */
  void print(ostream& os) const
  {
    os << "HEADERS,ELAPSED_MS" << endl;
    os << "VALUES," << getElapsedMS() << endl;
  }

  /**
   * Set the contents of the statistics table.  The idea is to keep one
   * table for each segment of code. New runs are appended to the
   * table.  
   * Returns: the new_row number if Ok 
   *          the setCell error code (< 0) if Error.
  */
  int setCellValues(CSVStringTable& st, unsigned row) const 
  {
    if (st.nRows() <= row) return -1;
    return st.setCell(row, "ELAPSED", getElapsedMS(), true);
  }

  void clearStats()
  {
    memset(&lap, 0, sizeof(lap));
    memset(&total, 0, sizeof(total));
  }

  double getElapsedMS() const
  {
    double elapsed;
    elapsed  = (total.tv_sec)  * 1000.0;   // sec to ms
    elapsed += (total.tv_usec) / 1000.0;   // us to ms
    return elapsed;
  }

 protected:

  timeval lap;
  timeval total;
};

/**
 * Class to record execution statistics. For performance reasons it
 * should be declared globally in order to ensure it is only
 * constructed and destructed once.
 */
class RunStatsRecorder
{
 public:

  /** Constructor. Create the RunStat items that should be measured. */
  RunStatsRecorder() : n_laps(0)
    {
#ifdef HAS_GETRUSAGE
      /* Add the resource usage statistics. */
      stat_items.push_back(new RUsageRunStat());
#endif
      /* Add the elapsed time statistic. */
      stat_items.push_back(new ElapsedTimeRunStat());
    }
  
  /** Destructor. */
  ~RunStatsRecorder()
    {
      vector<RunStat*>::iterator it;
      for (it=stat_items.begin(); it!=stat_items.end(); it++) {
	RunStat* i = *it;
	delete i;
      }
    }

  /** Starts recording the execution statistics. */
  void start()
  {
    vector<RunStat*>::iterator it;
    for (it=stat_items.begin(); it!=stat_items.end(); it++)
      (*it)->start();
  }

  /** Stops recording the execution statistics. */
  void stop()
  {
    vector<RunStat*>::reverse_iterator rit;
    for (rit=stat_items.rbegin(); rit!=stat_items.rend(); rit++)
      (*rit)->stop();
    n_laps++;
  }

  /** Stops the current measurement and starts a new one. */
  void lap() { stop(); start(); }

  /** Prints the statistics */
  void print(ostream& os) const
  {
    vector<RunStat*>::const_iterator it;
    for (it=stat_items.begin(); it!=stat_items.end(); it++)
      (*it)->print(os);
  }

  /**
   * Append metrics to the statistics table.  The idea is to keep one
   * table for each segment of code. New runs are appended to the
   * table.  
   * Returns: the new_row number if Ok 
   *          the setCell error code (< 0) if Error.
  */
  int append_to(CSVStringTable& st) const 
  {
    int ret;
    unsigned new_row = st.addRows(1);
    
    if ((ret=update_row(st, new_row)))
      return ret;
    else
      return new_row;
  }

  /**
   * Update a row at the statistics table.
   * Returns: 0 if Ok 
   *          the setCell error code (< 0) if Error.
  */
  int update_row(CSVStringTable& st, unsigned row) const 
  {
    int ret;

    vector<RunStat*>::const_iterator it;
    for (it=stat_items.begin(); it!=stat_items.end(); it++) {
      if ((ret=(*it)->setCellValues(st,row))) {
	return ret;
      }
    }
    
    return 0;
  }

  void clear() {
    
    vector<RunStat*>::iterator it;
    for (it=stat_items.begin(); it!=stat_items.end(); it++) 
      (*it)->clearStats();
    n_laps = 0;
  }
  
 private:

  /** Array of statistics items. */ 
  vector<RunStat*> stat_items;

  unsigned n_laps;  //N_LAPS

};

#endif // STATS_RECORDER

