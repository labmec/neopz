/**
 * @file
 * @brief Contains macros and functions to support execution statistics recording.
 */
#ifndef STATS_RECORDER_H

#ifndef VC
#define HAS_GETRUSAGE
#endif

#include <vector>
#include "csvtable.h"
#include "pzreal.h"
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

class PAPIRunStat : public RunStat
{
public:
    
    PAPIRunStat();
    
    /* Start recording the execution statistics. */
    void start() override;
    
    /* Stop recording the execution statistics. */
    void stop() override;
    
    /**
     * Print the statistics.
     * TODO: create a table, update it, and print the table.
     */
    void print(std::ostream& os) const override;
    
    /**
     * Append metrics to the statistics table.  The idea is to keep one
     * table for each segment of code. New runs are appended to the
     * table.
     * Returns: the new_row number if Ok
     *          the setCell error code (< 0) if Error.
     */
    int setCellValues(CSVStringTable& st, unsigned row) const override;
    
    void clearStats() override;
    
protected:
    
    float rtimeS, ptimeS, mflopsS;
    int64_t flpopsS;
    
    float rtimeACC, ptimeACC;
    int64_t flpopsACC;
    
};

    
#include<pz_gettime.h>
    
    class ElapsedTimeRunStat : public RunStat
    {
    public:
        
        ElapsedTimeRunStat();
        
        /* Start recording the execution statistics. */
        void start() override;
        
        /* Stop recording the execution statistics. */
        void stop() override;
        
        /**
         * Print the statistics
         */
        void print(std::ostream& os) const override;
        
        /**
         * Set the contents of the statistics table.  The idea is to keep one
         * table for each segment of code. New runs are appended to the
         * table.
         * Returns: the new_row number if Ok
         *          the setCell error code (< 0) if Error.
         */
        int setCellValues(CSVStringTable& st, unsigned row) const override;
        
        void clearStats() override;
        
        double getElapsedMS() const;
        
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
        RunStatsRecorder();
        
        /** Destructor. */
        ~RunStatsRecorder();
        
        /** Starts recording the execution statistics. */
        void start();
        
        /** Stops recording the execution statistics. */
        void stop();
        
        /** Stops the current measurement and starts a new one. */
        void lap();
        
        /** Prints the statistics */
        void print(std::ostream& os) const;
        
        /**
         * Append metrics to the statistics table.  The idea is to keep one
         * table for each segment of code. New runs are appended to the
         * table.
         * Returns: the new_row number if Ok
         *          the setCell error code (< 0) if Error.
         */
        int append_to(CSVStringTable& st) const;
        
        /**
         * Update a row at the statistics table.
         * Returns: 0 if Ok
         *          the setCell error code (< 0) if Error.
         */
        int update_row(CSVStringTable& st, unsigned row) const;
        
        void clear();
        
    private:
        
        /** Array of statistics items. */ 
        std::vector<RunStat*> stat_items;
        
        unsigned n_laps;  //N_LAPS
        
    };
    
#endif // STATS_RECORDER
    
