/**
 * @file
 * @brief Tests for decompose_ldlt
 * @author Edson Borin
 * @since 2012
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <cstdlib>

#include "pzbfilestream.h" // TPZBFileStream, TPZFileStream
#include "pzmd5stream.h"

#include <fstream>
#include <string>

#ifdef LOG4CXX
static LoggerPtr loggerconverge(Logger::getLogger("pz.converge"));
static LoggerPtr logger(Logger::getLogger("main"));
#endif

#include "pzskylmat.h"

//#include "timing_analysis.h"
#include "arglib.h"
#include "run_stats_table.h"

#ifdef HAS_GETRUSAGE
#include <sys/resource.h> // getrusage
#endif

#ifdef USING_TBB
#include "tbb/task_scheduler_init.h"
using namespace tbb; 
// If you have issues with: dyld: Library not loaded: libtbb.dylib
// try setting the LD path. Ex: 
//   export DYLD_FALLBACK_LIBRARY_PATH=/Users/borin/Desktop/neopz/tbb40_297oss/lib/
#endif

using namespace std;

void help(const char* prg)
{
    cout << "Compute the Decompose_LDLt method for the matrix" << endl;
    cout << endl;
    cout << "Usage: " << prg << "-if file [-v verbose_level] [-b] "
	 << "[-tot_rdt rdt_file] [-h]" << endl << endl;
    
    clarg::arguments_descriptions(cout, "  ", "\n");
} 

clarg::argString ifn("-ifn", "input matrix file name", "matrix.txt");
clarg::argInt verb_level("-v", "verbosity level", 0);
clarg::argBool bi("-b", "binary input file", false);
clarg::argBool h("-h", "help message", false);

/* Run statistics. */
RunStatsTable total_rst("-tot_rdt", 
			"Whole program (total) statistics raw data table");

class FileStreamWrapper
{
public: 
  FileStreamWrapper() {}
  ~FileStreamWrapper() {}
  
  void OpenWrite(const std::string& fn)
  {
    if (bi.was_set())
      bfs.OpenWrite(fn);
    else
      fs.OpenWrite(fn);
  }

  void OpenRead(const std::string& fn)
  {
    if (bi.was_set())
      bfs.OpenRead(fn);
    else
      fs.OpenRead(fn);
  }

  operator TPZStream&() 
  {
    if (bi.was_set())
      return bfs;
    else
      return fs;
  }

protected:

  TPZFileStream  fs;
  TPZBFileStream bfs;  
};


int main(int argc, char *argv[])
{
#ifdef USING_TBB
    task_scheduler_init init;
#endif
    
    /* Parse the arguments */
    if (clarg::parse_arguments(argc, argv)) {
        cerr << "Error when parsing the arguments!" << endl;
        return 1;
    }
    
    if (h.get_value() == true) {
        help(argv[0]);
        return 1;
    }
    
    /* Verbose macro. */
    unsigned verbose = verb_level.get_value();
#define VERBOSE(level,...) if (level <= verbose) cout << __VA_ARGS__
    
    if (verbose >= 1) {
        std::cout << "- Arguments -----------------------" << std::endl;
        clarg::values(std::cout, false);
        std::cout << "-----------------------------------" << std::endl;
    }

    /* Read the matrix. */
    TPZSkylMatrix<REAL> matrix;

    VERBOSE(1,"Reading input file: " << ifn.get_value() << std::endl);
    FileStreamWrapper input_file;
    input_file.OpenRead(ifn.get_value());
    matrix.Read(input_file,0);
    VERBOSE(1,"Reading input file: " << ifn.get_value() 
	    << " [DONE]" << std::endl);

    VERBOSE(1,"Starting Decompose_LDLt (Matrix dim = " << matrix.Dim() 
	    << ")" << std::endl);
    total_rst.start();
    matrix.Decompose_LDLt();
    total_rst.stop();
    VERBOSE(1,"Starting Decompose_LDLt [DONE]" << std::endl);

    // Return OK.
    return 0;
}
