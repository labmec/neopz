/**
 * @file
 * @brief  Tests for skyline matrices
 * @author Edson Borin
 * @since  2013
 */

#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include <iostream>
#include <cstdlib>

#include "pzbfilestream.h" // TPZBFileStream, TPZFileStream
#include "pzmd5stream.h"

#include <fstream>
#include <string>

#ifdef USING_TBB
#include "tbb/task_scheduler_init.h"
using namespace tbb;
#endif

#include "pzskylmat.h"

#include "arglib.h"
#include "run_stats_table.h"

using namespace std;

void help(const char* prg)
{
    cout << "Execute performance tests for skylmatrices" << endl;
    cout << endl;
    cout << "Usage: " << prg << "-op operation -[m|bm] matrixfn [-[m|bm]2 matrixfn] [-v level] [-perf_rdt rdt_file] [-h]";
    cout << "[extra_arguments] [-h]" << endl << endl;
    cout << "operation:" << endl;
    cout << " 0: dump skyline matrix statistics" << endl;
    cout << " 1: decompose matrix using Decompose_Cholesky()" << endl;
    cout << " 2: decompose matrix using Decompose_LDLt()" << endl;
    clarg::arguments_descriptions(cout, "  ", "\n");
}

clarg::argString m("-m", "input matrix file name (text format)", "matrix.txt");
clarg::argString bm("-bm", "input matrix file name (binary format)", "matrix.bin");
clarg::argString m2("-m2", "argument matrix file name (text format)", "matrix2.txt");
clarg::argString bm2("-bm2", "argument matrix file name (binary format)", "matrix2.bin");

clarg::argInt verb_level("-v", "verbosity level", 0);
int verbose = 0;
/* Verbose macro. */
#define VERBOSE(level,...) if (level <= verbose) cout << __VA_ARGS__

clarg::argInt mop("-op", "Matrix operation", 1);
clarg::argBool h("-h", "help message", false);
clarg::argString res_chk_t("-res_chk_t", "check the results using a reference data (text format)","");
clarg::argString res_chk_b("-res_chk_b", "check the results using a reference data (binary format)","");
clarg::argDouble res_chk_tol("-res_chk_tol", "error tolerance when checking results.", 1.e-12);
clarg::argString res_dump_t("-res_dump_t", "write final results to a text file.", "dump_matrix.txt");
clarg::argString res_dump_b("-res_dump_b", "write final results to a binary file.", "dump_matrix.bin");

/* Run statistics. */
RunStatsTable operation_rst("-perf_rdt", "Raw data table file to add matrix operation performance statistics");

int run_decompose_cholesky();
int run_decompose_ldlt();
int dump_matrix_stats();

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
    
    verbose = verb_level.get_value();
    if (verbose >= 1) {
        cout << "- Arguments -----------------------" << endl;
        clarg::values(cout, false);
        cout << "-----------------------------------" << endl;
    }
    
    if (h.get_value() == true) {
        help(argv[0]);
        return 1;
    }
    
    if (!mop.was_set()) {
        cerr << "You must provide the operation -op." << endl;
        help(argv[0]);
        return 1;
    }
    
    int status = 0; // Ok.
    
    switch (mop.get_value()) {
        case 0:
            status = dump_matrix_stats();
            break;
        case 1: 
            status = run_decompose_cholesky(); 
            break;
        case 2: 
            status = run_decompose_ldlt(); 
            break;
        default:
            cerr << "ERROR: Invalid matrix operation type." << endl;
            status = 1;
    }
    
    if (status != 0) {
        cerr << "ERROR when executing the experiment." << endl;
    }
    
    return status;
}

// -- Helper functions ---------------------------------------

class FileStreamWrapper
{
public:
    FileStreamWrapper(bool b) : binary(b)
    {}
    ~FileStreamWrapper() {}
    
    void OpenWrite(const string& fn)
    {
        if (binary)
            bfs.OpenWrite(fn);
        else
            fs.OpenWrite(fn);
    }
    
    void OpenRead(const string& fn)
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

int res_check(TPZSkylMatrix<REAL>& matrix) 
{
    string filename = 0;
    bool binary = false;
    if (res_chk_t.was_set()) {
        if (res_chk_b.was_set()) {
            cerr << "Both -res_chk_t and -res_chk_b were set. Use only one option." << endl;
            return 1;
        }
        else {
            filename = res_chk_t.get_value();
        }
    }
    else {
        binary = true;
        filename = res_chk_b.get_value();
    }
    
    VERBOSE(1,"Checking result using reference matrix : " << filename << endl);
    FileStreamWrapper file(binary);
    file.OpenRead(filename);
    /* Reference matrix. */
    TPZSkylMatrix<REAL> ref_matrix;
    ref_matrix.Read(file,0);
    
    int max_j = matrix.Cols();
    if (max_j != ref_matrix.Cols()) {
        cerr << "Result matrix has " << max_j
        << " cols while reference matrix has "
        << ref_matrix.Cols() << endl;
        VERBOSE(1,"Checking result using reference matrix : " << filename << "[FAILED]" << endl);
        return 1;
    }
    
    REAL error_tolerance = res_chk_tol.get_value();
    REAL max_error = 0.0;
    int ret = 0; // OK
    for (int j=0; j<max_j; j++) {
        int col_height = matrix.SkyHeight(j);
        if (col_height != ref_matrix.SkyHeight(j)) {
            cerr << "Column " << j << " of result matrix has " << col_height
            << " non zero rows while reference matrix has "
            << ref_matrix.SkyHeight(j) << endl;
            VERBOSE(1,"Checking result using reference matrix : " << filename << "[FAILED]" << endl);
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
        VERBOSE(1,"Checking result using reference matrix : " << filename << "[FAILED]" << endl);
    }
    else {
        VERBOSE(1,"Checking result using reference matrix : " << filename << "[OK]" << endl);
    }
    return ret;
}

int res_dump(TPZSkylMatrix<REAL>& matrix) 
{
    if (!res_dump_t.was_set() && !res_dump_b.was_set())
        return 0; // nothing to dump
    
    string filename = 0;
    bool binary = false;
    if (res_dump_t.was_set()) {
        if (res_dump_b.was_set()) {
            cerr << "Both -res_dump_t and -res_dump_b were set. Use only one option." << endl;
            return 1;
        }
        else {
            filename = res_dump_t.get_value();
        }
    }
    else {
        binary = true;
        filename = res_dump_b.get_value();
    }
    
    VERBOSE(1,"Dumping result to : " << filename << endl);
    FileStreamWrapper file(binary);
    file.OpenWrite(filename);
    matrix.Write(file,0);
    VERBOSE(1,"Dumping result to : " << filename << "[DONE]" << endl);
    return 0;
}

int read_input_matrix(TPZSkylMatrix<REAL>& matrix) 
{
    string inputfn = 0;
    bool binary = false;
    if (!m.was_set() && !bm.was_set()) {
        cerr << "Please provide an input matrix using -m or -bm." << endl;
        return 1;
    }
    if (m.was_set()) {
        if (bm.was_set()) {
            cerr << "Both -m and -bm were set. Use only one option." << endl;
            return 1;
        }
        else {
            inputfn = m.get_value();
        }
    }
    else {
        binary = true;
        inputfn = bm.get_value();
    }
    
    VERBOSE(1,"Reading input file: " << inputfn << endl);
    FileStreamWrapper input_file(binary);
    input_file.OpenRead(inputfn);
    matrix.Read(input_file,0);
    VERBOSE(1,"Reading input file: " << inputfn << "[DONE]" << endl);
    return 0;
}

// -- dump_matrix_stats ---------------------------------------
int dump_matrix_stats() 
{
    TPZSkylMatrix<REAL> matrix;
    
    if (read_input_matrix(matrix)) 
        return 1;
    
    unsigned n = matrix.Dim();
    uint64_t n_sky_items = 0;
    uint64_t max_height = 0;
    for (unsigned i=0; i<n; i++) {
        unsigned height = matrix.SkyHeight(i);
        if (verbose > 2) {
            cout << "col " << i << " height = " << height << endl;
        }
        n_sky_items += height;
        if (height > max_height) max_height = height;
    }
    uint64_t n2 = n * n;
    double av_height = (double) n_sky_items / (double) n;
    cout << "N         = " << n << endl;
    cout << "N^2       = " << n2 << endl;
    cout << "Sky items = " << n_sky_items << endl;
    cout << "N^2 / Sky items = " << (double) n2 / (double) n_sky_items << endl;
    cout << "Avg. Height = " << av_height << endl;
    cout << "Max. Height = " << max_height << endl;
    
    return 0;
} 


// -- run_decompose_cholesky ---------------------------------------
int run_decompose_cholesky() 
{
    TPZSkylMatrix<REAL> matrix;
    
    if (read_input_matrix(matrix)) 
        return 1;
    
    operation_rst.start();
    //  int ret = matrix.Decompose_Cholesky();
    operation_rst.stop();
    
    if (res_dump(matrix)) 
        return 1;
    
    if (res_chk_t.was_set() || res_chk_b.was_set())
        return res_check(matrix);
    else
        return 0;
}

// -- run_decompose_ldlt ---------------------------------------
int run_decompose_ldlt() 
{
    TPZSkylMatrix<REAL> matrix;
    
    if (read_input_matrix(matrix)) 
        return 1;
    
    operation_rst.start();
    //    int ret = matrix.Decompose_LDLt();
    operation_rst.stop();
    
    if (res_dump(matrix)) 
        return 1;
    
    if (res_chk_t.was_set() || res_chk_b.was_set())
        return res_check(matrix);
    else
        return 0;
}
