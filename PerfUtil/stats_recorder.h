/**
 * @file
 * @brief Contains macros and functions to support execution statistics recording.
 */
#ifndef STATS_RECORDER_H

#ifndef VC
#define HAS_GETRUSAGE
#endif

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
 
 / * Contains statistics
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
    
    virtual ~RunStat()
    {
        
    }
    /**
     * Print the statistics
     */
    virtual void print(std::ostream& os) const = 0;
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

#ifdef USING_PAPI
#include <papi.h>
class PAPIRunStat : public RunStat
{
public:
    
    PAPIRunStat()
    {
        clearStats();
    }
    
    /* Start recording the execution statistics. */
    void start()
    {
        PAPI_flops ( &rtimeS, &ptimeS, &flpopsS, &mflopsS);
    }
    
    /* Stop recording the execution statistics. */
    void stop()
    {
        float rtimeP, ptimeP, mflopsP;
        int64_t flpopsP;
        PAPI_flops ( &rtimeP, &ptimeP, &flpopsP, &mflopsP);
        rtimeACC += (rtimeP-rtimeS);
        ptimeACC += (ptimeP-ptimeS);
        flpopsACC += (flpopsP-flpopsS);
    }
    
    /**
     * Print the statistics.
     * TODO: create a table, update it, and print the table.
     */
    void print(std::ostream& os) const
    {
        os << "HEADERS,PAPI_RTIME,PAPI_PTIME,PAPI_FLPOPS,PAPI_MFLOPS" << endl;
        os << "VALUES,"
        << rtimeACC << ","
        << ptimeACC << ","
        << flpopsACC << ","
        << ((double) flpopsACC  / (double) ptimeACC)/1000000.0 << endl;
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
        int ret;
        if ((ret = st.setCell(row,"PAPI_RTIME",rtimeACC,true)))
            return ret;
        if ((ret = st.setCell(row,"PAPI_PTIME",ptimeACC,true)))
            return ret;
        if ((ret = st.setCell(row,"PAPI_FLPOPS",flpopsACC,true)))
            return ret;
        double av_mflops = ((double) flpopsACC / (double) ptimeACC) / 1000000.0;
        if ((ret = st.setCell(row,"PAPI_AV_MFLOPS",av_mflops,true)))
            return ret;
        
        return 0; // Return ok
    }
    
    void clearStats()
    {
        memset(&rtimeACC, 0, sizeof(rtimeACC));
        memset(&ptimeACC, 0, sizeof(ptimeACC));
        memset(&flpopsACC, 0, sizeof(flpopsACC));
    }
    
protected:
    
    float rtimeS, ptimeS, mflopsS;
    int64_t flpopsS;
    
    float rtimeACC, ptimeACC;
    int64_t flpopsACC;
    
};

#endif

#ifdef REALpzfpcounter
#include "pzreal.h"
class PZFPCountRunStat : public RunStat
{
public:
    
    PZFPCountRunStat()
    {
        clearStats();
    }
    
    /* Start recording the execution statistics. */
    void start()
    {
        start_counter.copy(TPZFlopCounter::gCount);
    }
    
    /* Stop recording the execution statistics. */
    void stop()
    {
        TPZCounter stop_counter;
        stop_counter.copy(TPZFlopCounter::gCount);
        stop_counter -= start_counter;
        acc_counter += stop_counter;
    }
    
    /**
     * Print the statistics.
     * TODO: create a table, update it, and print the table.
     */
    void print(ostream& os) const
    {
        stringstream headers;
        stringstream values;
        
        headers << "HEADERS";
        for (int i=0; i<gNumOp; i++) {
            headers << "," << "PZCOUNT_" << OpNames[i];
        }
        values << "VALUES";
        for (int i=0; i<gNumOp; i++) {
            values << "," << acc_counter[i];
        }
        os << headers << endl;
        os << values << endl;
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
        int ret;
        
        for (int i=0; i<gNumOp; i++) {
            stringstream header;
            header << "PZCOUNT_" << OpNames[i];
            if ((ret = st.setCell(row,header,acc_counter[i],true)))
                return ret;
        }
        
        return 0; // Return ok
    }
    
    void clearStats()
    {
        acc_counter.clear();
    }
    
protected:
    
    TPZCounter start_counter;
    TPZCounter acc_counter;
    
};

#endif

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
        
#define SET_TOTAL_TIMEVAL(fld)						\
total_self.fld.tv_sec += (self.fld.tv_sec - lap_self.fld.tv_sec); \
total_self.fld.tv_usec += (self.fld.tv_usec - lap_self.fld.tv_usec); \
total_children.fld.tv_sec += (children.fld.tv_sec - lap_children.fld.tv_sec); \
total_children.fld.tv_usec += (children.fld.tv_usec - lap_children.fld.tv_usec)
        
        SET_TOTAL_TIMEVAL(ru_utime);   /* user time used */
        SET_TOTAL_TIMEVAL(ru_stime);   /* system time used */
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
    void print(std::ostream& os) const
    {
        std::stringstream header;
        std::stringstream values;
        
#define PRINT_FLD(hd,fld) header << ",SELF_" << hd << ",CHD_" << hd; values << "," << total_self.fld << "," << total_children.fld
        
#define TIMEVAL_TO_DOUBLE_MS(s) ( ((double) s.tv_sec) * 1000.0 + ((double) s.tv_usec) / 1000.0)
        
#define PRINT_TIMEVAL_FLD(hd,fld)					\
header << ",SELF_" << hd << ",CHD_" << hd;			\
values << "," << TIMEVAL_TO_DOUBLE_MS(total_self.fld) << "," << TIMEVAL_TO_DOUBLE_MS(total_children.fld)
        
        PRINT_TIMEVAL_FLD("RU_UTIME",ru_utime);           /* user time used */
        PRINT_TIMEVAL_FLD("RU_STIME",ru_stime);           /* system time used */
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
        os << "HEADERS" << header.str() << std::endl;
        os << "VALUES" << values.str() << std::endl;
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
{ int ret; std::string header("SELF_"); header += hd;			\
if ((ret = st.setCell(row, header, total_self.fld, true)))	\
return ret;							\
header = "CHD_"; header += hd;					\
if ((ret = st.setCell(row, header, total_children.fld, true)))	\
return ret;							\
}
        
#define APPEND_TIMEVAL_FLD(hd,fld)						\
{ int ret; std::string header("SELF_"); header += hd;			\
if ((ret = st.setCell(row, header, TIMEVAL_TO_DOUBLE_MS(total_self.fld), true))) \
return ret;							\
header = "CHD_"; header += hd;					\
if ((ret = st.setCell(row, header, TIMEVAL_TO_DOUBLE_MS(total_children.fld), true)))	\
return ret;							\
}
        APPEND_TIMEVAL_FLD("RU_UTIME",ru_utime);           /* user time used */
        APPEND_TIMEVAL_FLD("RU_STIME",ru_stime);           /* system time used */
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

#ifdef HAS_LIBPFM
//#include <perfmon/pfmlib_perf_event.h>
//#include "perf_util.h"

/*
 TODO:
 1- Initialize the libpfm
 ret = pfm_initialize();
 if (ret != PFM_SUCCESS) error;
 
 */

class LIBPFMCounter
{
Public:
    
private:
    perf_event_desc_t event;
    
};

class LIBPFMRunStat : public RunStat
{
private:
    
    bool counters_enabled;
    bool counters_created;
    vector<perf_event_desc_t> counters;
    
public:
    
    LIBPFMRunStat() :
    counters_enable(false), counters_created(false)
    {
        clearStats();
    }
    
    /* Start recording the execution statistics. */
    void start()
    {
        // Check whether libpfm was properly initialized.
        if (!libpfm_man.ok()) return;
        
        // Create counters
        for (int i=0; i<counters.size(); i++) {
            fds[i].hw.read_format = PERF_FORMAT_SCALE;
            fds[i].hw.disabled = 1; /* do not start now */
            /* each event is in an independent group (multiplexing likely) */
            fds[i].fd = perf_event_open(&fds[i].hw, 0, -1, -1, 0);
            if (fds[i].fd == -1)
                err(1, "cannot open event %d", i);
        }
        
        counters_created;
        
        // Enable counters attached to this thread
        if (prctl(PR_TASK_PERF_EVENTS_ENABLE)) {
            // ERROR(prctl(enable) failed)
            
        }
        counters_enable = true;
    }
    
    /* Stop recording the execution statistics. */
    void stop()
    {
        if (counters_enable != true) return;
        if (prctl(PR_TASK_PERF_EVENTS_DISABLE)) {
            // ERROR (prctl(disable) failed)
            
            
            
            //struct rusage self, children;
            //getrusage(RUSAGE_SELF, &self);
            //getrusage(RUSAGE_CHILDREN, &children);
            //#define SET_TOTAL(fld) total_self.fld += (self.fld - lap_self.fld); total_children.fld += (children.fld - lap_children.fld)
            ////struct timeval ru_utime; /* user time used */
            ////struct timeval ru_stime; /* system time used */
            //SET_TOTAL(ru_maxrss);          /* integral max resident set size */
            //SET_TOTAL(ru_ixrss);           /* integral shared text memory size */
            //SET_TOTAL(ru_idrss);           /* integral unshared data size */
            //SET_TOTAL(ru_isrss);           /* integral unshared stack size */
            //SET_TOTAL(ru_minflt);          /* page reclaims */
            //SET_TOTAL(ru_majflt);          /* page faults */
            //SET_TOTAL(ru_nswap);           /* swaps */
            //SET_TOTAL(ru_inblock);         /* block input operations */
            //SET_TOTAL(ru_oublock);         /* block output operations */
            //SET_TOTAL(ru_msgsnd);          /* messages sent */
            //SET_TOTAL(ru_msgrcv);          /* messages received */
            //SET_TOTAL(ru_nsignals);        /* signals received */
            //SET_TOTAL(ru_nvcsw);           /* voluntary context switches */
            //SET_TOTAL(ru_nivcsw);          /* involuntary context switches */
        }
        
        /**
         * Print the statistics.
         * TODO: create a table, update it, and print the table.
         */
        void print(ostream& os) const
        {
            //    stringstream header;
            //    stringstream values;
            //#define PRINT_FLD(hd,fld) header << ",SELF_" << hd << ",CHD_" << hd; values << "," << total_self.fld << "," << total_children.fld
            //
            //    //struct timeval ru_utime; /* user time used */
            //    //struct timeval ru_stime; /* system time used */
            //    PRINT_FLD("RU_MAXRSS",ru_maxrss);          /* integral max resident set size */
            //    PRINT_FLD("RU_IXRSS",ru_ixrss);           /* integral shared text memory size */
            //    PRINT_FLD("RU_IDRSS",ru_idrss);           /* integral unshared data size */
            //    PRINT_FLD("RU_ISRSS",ru_isrss);           /* integral unshared stack size */
            //    PRINT_FLD("RU_MINFLT",ru_minflt);          /* page reclaims */
            //    PRINT_FLD("RU_MAJFLT",ru_majflt);          /* page faults */
            //    PRINT_FLD("RU_NSWAP",ru_nswap);           /* swaps */
            //    PRINT_FLD("RU_INBLOCK",ru_inblock);         /* block input operations */
            //    PRINT_FLD("RU_OUBLOCK",ru_oublock);         /* block output operations */
            //    PRINT_FLD("RU_MSGND",ru_msgsnd);          /* messages sent */
            //    PRINT_FLD("RU_MSGRCV",ru_msgrcv);          /* messages received */
            //    PRINT_FLD("RU_NSIGNAL",ru_nsignals);        /* signals received */
            //    PRINT_FLD("RU_NVCSW",ru_nvcsw);           /* voluntary context switches */
            //    PRINT_FLD("RU_NIVSW",ru_nivcsw);          /* involuntary context switches */
            //    os << "HEADERS" << header.str() << std::endl;
            //    os << "VALUES" << values.str() << std::endl;
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
            //    if (st.nRows() <= row) return -1;
            //
            //#define APPEND_LONG_FLD(hd,fld)						\
            //    { int ret; string header("SELF_"); header += hd;			\
            //      if ((ret = st.setCell(row, header, total_self.fld, true)))	\
            //	return ret;							\
            //      header = "CHD_"; header += hd;					\
            //      if ((ret = st.setCell(row, header, total_children.fld, true)))	\
            //	return ret;							\
            //    }
            //
            //    //struct timeval ru_utime; /* user time used */
            //    //struct timeval ru_stime; /* system time used */
            //    APPEND_LONG_FLD("RU_MAXRSS",ru_maxrss);          /* integral max resident set size */
            //    APPEND_LONG_FLD("RU_IXRSS",ru_ixrss);           /* integral shared text memory size */
            //    APPEND_LONG_FLD("RU_IDRSS",ru_idrss);           /* integral unshared data size */
            //    APPEND_LONG_FLD("RU_ISRSS",ru_isrss);           /* integral unshared stack size */
            //    APPEND_LONG_FLD("RU_MINFLT",ru_minflt);          /* page reclaims */
            //    APPEND_LONG_FLD("RU_MAJFLT",ru_majflt);          /* page faults */
            //    APPEND_LONG_FLD("RU_NSWAP",ru_nswap);           /* swaps */
            //    APPEND_LONG_FLD("RU_INBLOCK",ru_inblock);         /* block input operations */
            //    APPEND_LONG_FLD("RU_OUBLOCK",ru_oublock);         /* block output operations */
            //    APPEND_LONG_FLD("RU_MSGND",ru_msgsnd);          /* messages sent */
            //    APPEND_LONG_FLD("RU_MSGRCV",ru_msgrcv);          /* messages received */
            //    APPEND_LONG_FLD("RU_NSIGNAL",ru_nsignals);        /* signals received */
            //    APPEND_LONG_FLD("RU_NVCSW",ru_nvcsw);           /* voluntary context switches */
            //    APPEND_LONG_FLD("RU_NIVSW",ru_nivcsw);          /* involuntary context switches */
            //
            return 0; // Return ok
        }
        
        void clearStats()
        {
            //    memset(&lap_self, 0, sizeof(lap_self));
            //    memset(&lap_children, 0, sizeof(lap_children));
            //    memset(&total_self, 0, sizeof(total_self));
            //    memset(&total_children, 0, sizeof(total_children));
        }
        
    protected:
        
        //  struct rusage lap_self, lap_children;
        //  struct rusage total_self, total_children;
        
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
        void print(std::ostream& os) const
        {
            os << "HEADERS,ELAPSED_MS" << std::endl;
            os << "VALUES," << getElapsedMS() << std::endl;
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
#ifdef HAS_LIBPFM
            /* Add the resource usage statistics. */
            stat_items.push_back(new LIBPFMRunStat());
#endif
#ifdef HAS_GETRUSAGE
            /* Add the resource usage statistics. */
            stat_items.push_back(new RUsageRunStat());
#endif
#ifdef USING_PAPI
            /* Add the PAPI counters statistics. */
            stat_items.push_back(new PAPIRunStat());
#endif
            /* Add the elapsed time statistic. */
            stat_items.push_back(new ElapsedTimeRunStat());
        }
        
        /** Destructor. */
        ~RunStatsRecorder()
        {
            std::vector<RunStat*>::iterator it;
            for (it=stat_items.begin(); it!=stat_items.end(); it++) {
                RunStat* i = *it;
                delete i;
            }
            /** @brief vector::clear Removes all elements from the vector (which are destroyed), leaving the container with a size of 0. */
            stat_items.clear();
        }
        
        /** Starts recording the execution statistics. */
        void start()
        {
            std::vector<RunStat*>::iterator it;
            for (it=stat_items.begin(); it!=stat_items.end(); it++)
                (*it)->start();
        }
        
        /** Stops recording the execution statistics. */
        void stop()
        {
            std::vector<RunStat*>::reverse_iterator rit;
            for (rit=stat_items.rbegin(); rit!=stat_items.rend(); rit++)
                (*rit)->stop();
            n_laps++;
        }
        
        /** Stops the current measurement and starts a new one. */
        void lap() { stop(); start(); }
        
        /** Prints the statistics */
        void print(std::ostream& os) const
        {
            std::vector<RunStat*>::const_iterator it;
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
            
            std::vector<RunStat*>::const_iterator it;
            for (it=stat_items.begin(); it!=stat_items.end(); it++) {
                if ((ret=(*it)->setCellValues(st,row))) {
                    return ret;
                }
            }
            
            return 0;
        }
        
        void clear() {
            
            std::vector<RunStat*>::iterator it;
            for (it=stat_items.begin(); it!=stat_items.end(); it++)
                (*it)->clearStats();
            n_laps = 0;
        }
        
    private:
        
        /** Array of statistics items. */ 
        std::vector<RunStat*> stat_items;
        
        unsigned n_laps;  //N_LAPS
        
    };
    
#endif // STATS_RECORDER
    
