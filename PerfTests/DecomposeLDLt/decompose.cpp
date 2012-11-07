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
	 << "[-tot_rdt rdt_file] [-op matrix_operation] [-h]" << endl << endl;
    cout << "matrix_operation:" << endl;
    cout << " 0: Decompose_LDLt()" << endl;
    cout << " 1: Decompose_LDLt2()" << endl;
    cout << " 2: Decompose_Cholesky()" << endl;
    clarg::arguments_descriptions(cout, "  ", "\n");
} 

clarg::argString ifn("-ifn", "input matrix file name (use -bi to read from binary files)", "matrix.txt");
clarg::argInt verb_level("-v", "verbosity level", 0);
clarg::argInt mop("-op", "Matrix operation", 1);
clarg::argBool br("-br", "binary reference. Reference decomposed matrix file format == binary.", false);
clarg::argBool bi("-bi", "binary input. Input file format == binary.", false);
clarg::argBool bd("-bd", "binary dump. Dump file format == binary.", false);
clarg::argBool h("-h", "help message", false);
clarg::argInt mstats("-mstats verbosity", "Matrix statistics vebosity level.", 0);
clarg::argString gen_dm_sig("-gen_dm_md5", "generates MD5 signature for decomposed matrix into file.", "decomposed_matrix.md5");
clarg::argString chk_dm_sig("-chk_dm_md5", "compute MD5 signature for decomposed matrix and check against MD5 at file.", "decomposed_matrix.md5");
clarg::argString chk_dm_error("-chk_dm_error", "check the decomposed matrix error against a reference matrix. (use -br to read from binary files)", "ref_decomposed_matrix.txt");
clarg::argDouble error_tol("-error_tol", "error tolerance.", 1.e-12);
clarg::argString dump_dm("-dump_dm", "dump decomposed matrix. (use -bd for binary format)", "dump_matrix.txt");

/* Run statistics. */
RunStatsTable total_rst("-tot_rdt", 
			"Whole program (total) statistics raw data table");

class FileStreamWrapper
{
public: 
  FileStreamWrapper(bool b) : binary(b)
    {}
  ~FileStreamWrapper() {}
  
  void OpenWrite(const std::string& fn)
  {
    if (binary)
      bfs.OpenWrite(fn);
    else
      fs.OpenWrite(fn);
  }

  void OpenRead(const std::string& fn)
  {
    if (binary)
      bfs.OpenRead(fn);
    else
      fs.OpenRead(fn);
  }

  operator TPZStream&() 
  {
    if (binary)
      return bfs;
    else
      return fs;
  }

protected:

  bool binary;
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
    FileStreamWrapper input_file(bi.get_value());
    input_file.OpenRead(ifn.get_value());
    matrix.Read(input_file,0);
    VERBOSE(1,"Reading input file: " << ifn.get_value() 
	    << " [DONE]" << std::endl);

#define CASE_OP(opid,method)				\
    case opid:						\
      total_rst.start();				\
      matrix.method;					\
      total_rst.stop();					\
      break

    switch (mop.get_value()) {
      CASE_OP(0,Decompose_LDLt());
      CASE_OP(1,Decompose_LDLt2());
      CASE_OP(2,Decompose_Cholesky());
    default:
      std::cerr << "ERROR: Invalid matrix operation type." << std::endl;
    }

    if (mstats.get_value() > 0) {
      unsigned n = matrix.Dim();
      unsigned long long n_sky_items = 0;
      unsigned long long max_height = 0;
      for (unsigned i=0; i<n; i++) {
	unsigned height = matrix.SkyHeight(i);
	if (mstats.get_value() > 1) {
	  cout << "col " << i << " height = " << height << endl;
	}
	n_sky_items += height;
	if (height > max_height) max_height = height;
      }
      unsigned long long n2 = n * n;
      double av_height = (double) n_sky_items / (double) n;
      cout << "N         = " << n << endl;
      cout << "N^2       = " << n2 << endl;
      cout << "Sky items = " << n_sky_items << endl;
      cout << "N^2 / Sky items = " << (double) n2 / (double) n_sky_items << endl;
      cout << "Avg. Height = " << av_height << endl;
      cout << "Max. Height = " << max_height << endl;
    }

    /** Dump decomposed matrix */
    if (dump_dm.was_set()) {
      VERBOSE(1, "Dumping decomposed matrix into: " << 
	      dump_dm.get_value() << endl);
      FileStreamWrapper dump_file(bd.get_value());
      dump_file.OpenWrite(dump_dm.get_value());
      matrix.Write(dump_file, 0);
    }

    /* Gen/Check MD5 signature */
    if (gen_dm_sig.was_set() || chk_dm_sig.was_set()) {
      TPZMD5Stream sig;
      matrix.Write(sig, 1);
      int ret;
      if (chk_dm_sig.was_set()) {
	if ((ret=sig.CheckMD5(chk_dm_sig.get_value()))) {
	  cerr << "ERROR(ret=" << ret << ") : MD5 Signature for "
	       << "decomposed matrixdoes not match." << endl;
	  return 1;
	}
	else {
	  cout << "Checking decomposed matrix MD5 signature: [OK]" << endl;
	}
      }
      if (gen_dm_sig.was_set()) {
	if ((ret=sig.WriteMD5(gen_dm_sig.get_value()))) {
	  cerr << "ERROR (ret=" << ret << ") when writing the "
	       << "decomposed matrix MD5 signature to file: " 
	       << gen_dm_sig.get_value() << endl;
	  return 1;
	}
      }
    }

    int ret=0; // Ok

    /** Check decomposed matrix */
    if (chk_dm_error.was_set()) {
      VERBOSE(1, "Checking decomposed matrix error: " << 
	      chk_dm_error.get_value() << endl);
      FileStreamWrapper ref_file(br.get_value());
      ref_file.OpenRead(chk_dm_error.get_value());
      /* Reference matrix. */
      TPZSkylMatrix<REAL> ref_matrix;
      ref_matrix.Read(ref_file,0);
      int max_j = matrix.Cols();
      if (max_j != ref_matrix.Cols()) {
	cerr << "Decomposed matrix has " << max_j
	     << " cols while reference matrix has "
	     << ref_matrix.Cols() << endl;
	return 1;
      }
      REAL error_tolerance = error_tol.get_value();
      REAL max_error = 0.0;
      for (int j=0; j<max_j; j++) {
	int col_height = matrix.SkyHeight(j);
	if (col_height != ref_matrix.SkyHeight(j)) {
	  cerr << "Column " << j << " of decomposed matrix has " << col_height
	       << " non zero rows while reference matrix has "
	       << ref_matrix.SkyHeight(j) << endl;
	  return 1;
	}
	int min_i = (j+1) - col_height;
	for (int i=min_i; i<=j; i++) {

	  REAL dm_ij = matrix.s(i,j);
	  REAL rm_ij = ref_matrix.s(i,j);
	  if (dm_ij != rm_ij) {
	    REAL diff = abs(dm_ij - rm_ij);
	    if (diff >= error_tolerance) {
	      VERBOSE(1, "diff(" << diff << ") tolerance (" << error_tolerance 
		      << "). dm[" << i << "][" << j << "] (" << dm_ij
		      << ") != rm[" << i << "][" << j << "] (" << rm_ij 
		      << ")." << endl);
	      ret = 1;
	      max_error = (max_error < diff)?diff:max_error;
	    }
	  }
	}
      }
      if (ret != 0) {
	cerr << "Error ("<< max_error <<") > error tolerance ("
	     << error_tolerance <<")" <<  endl;
      }
  }

    return ret;
}
