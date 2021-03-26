/**
 * @file
 * @brief Contains a class to record running statistics on CSV tables.
 */
#ifndef RUN_STATS_TABLE_H

#include"stats_recorder.h"
#include <fstream>

#include"arglib.h"

/**
 * Class to record execution statistics into tables. For performance
 * reasons it should be declared globally in order to ensure it is
 * only constructed and destructed once.
 *
 * Typical usage
 *
 *   RunStatsTable dohrstruct_create_rst("-cre_rst", "Run statistics table for dohrstruct->Create()"za);
 *
 *   // Track statistics. Only measure if "-cre_rst" was set.
 *   dohrstruct_create_rst.start();
 *   dohrstruct_create_rst.lap();
 *   dohrstruct_create_rst.stop();
 *
 *   // Optional: Updates the table contents. Reads the file into the
 *   // table, add the statistics and writes the file back. Only performed
 *   // if the filename was provided.
 *   dohrstruct_create_rst.flush_to_file();
 *
 *   // Automatically calls flush_to_file().
 *   ~dohrstruct_create_srt()
 */
class RunStatsTable
{
public:
    
    /** Constructor. */
    RunStatsTable(const char* arg, const char* desc) : filename(arg,desc,""), stats_recorded(false)
    {}
    
    /** Destructor. Update the file. */
    ~RunStatsTable()
    {
        if (filename.was_set()) {
            switch(flush_to_file())
            {
                case 1: std::cerr << "WARNING: number of cells in one or more rows does not "
                    << "match the number of column headers when reading the \""
                    << filename.get_value() << "\" run stats file." << std::endl; break;
                case -1: std::cerr << "WARNING: could not read the table headers when reading the \""
                    << filename.get_value() << "\" run stats file."; break;
                case -2: std::cerr << "WARNING: could not open the \"" << filename.get_value()
                    << "\" run stats file for read." << std::endl; break;
                case -3: std::cerr << "WARNING: could not open the \"" << filename.get_value()
                    << "\" run stats file for write." << std::endl; break;
                case -4: std::cerr << "WARNING: error when appending statistics to \"" << filename.get_value()
                    << "\" run stats file for write." << std::endl; break;
                default:
                    break;
            }
        }
    }
    
    /** Track statistics. Starts the measurement. */
    void start() {if (filename.was_set()) { stats_recorded=true; stats.start(); }}
    /** Track statistics. Add lap. */
    void lap()   {if (filename.was_set()) { stats_recorded=true; stats.lap(); }}
    /** Track statistics. Stops the measurement. */
    void stop()  {if (filename.was_set()) { stats_recorded=true; stats.stop(); }}
    
    /** Optional: Updates the table contents. Reads the file into the
     *  table, add the statistics and writes the file back. Only performed
     *  if the filename was provided.
     *  Returns  0 if ok
     *           1 if number of cells in one or more rows does not match
     *             the number of column headers when reading the file.
     *           2 nothing to be flushed.
     *          -1 if could not read the table headers
     *          -2 could not open the file for read.
     *          -3 could not open the file for write.
     *          -4 error when appending statistics to the table.
     **/
    int flush_to_file()
    {
        CSVStringTable table;
        int ret;
        
        if (!stats_recorded) return 2;
        
        /* Read the table. */
        if ((ret=read_from_file(table, filename.get_value()))) {
            if (ret != -2) return ret;
            // Continue if could not open for read (-2).
        }
        
        /* Update the table. */
        if (stats.append_to(table) < 0) {
            return -4;
        }
        
        stats.clear();
        stats_recorded = false;
        
        /* Write the table back. */
        return write_to_file(table, filename.get_value());
    }
    
    bool was_set() const {
        return filename.was_set();
    }
    
private:
    
    /** Reads the table contents from a file.
     *  Returns  0 if ok,
     *           1 if number of cells in one or more rows does not match the number of column headers
     *          -1 if could not read the table headers
     *          -2 could not open the file.
     */
    int read_from_file(CSVStringTable& table, std::string filename)
    {
        std::ifstream ifs(filename.c_str());
        
        if (ifs.is_open()) {
            int ret=table.read(ifs);
            ifs.close();
            return ret;
        }
        else
            return -2;
    }
    
    /** Writes the table contents to a file.
     *  Returns  0 if ok.
     *          -3 could not open the file for write.
     */
    int write_to_file(const CSVStringTable& table, std::string filename)
    {
        std::ofstream ofs(filename.c_str());
        
        if (ofs.is_open()) {
            table.write(ofs);
            ofs.close();
            return 0;
        }
        else {
            return -3;
        }
    }
    
    RunStatsRecorder stats;
    clarg::argString filename;
    bool stats_recorded;
};


#endif // RUN_STATS_TABLE

