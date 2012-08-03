/**
 * @file
 * @brief Measure the amount of GFlops/s the system can sustain.
 * @author Edson Borin
 * @since 2012
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>

//#include "timing_analysis.h"
#include "arglib.h"
#include "run_stats_table.h"

using namespace std;

void help(const char* prg)
{
    cout << "Perform an array multiplication: a[1...n] = b[1...n] x c[1...n]"
         << endl;
    
    cout << "Usage: " << prg << "[-asz N] [-modb N] [-modc N] [-moda N] [-nt N]" 
         << endl << endl;
    
    clarg::arguments_descriptions(cout, "  ", "\n");
} 

clarg::argInt  asz("-asz", "array size", 10000000);
clarg::argInt  moda("-moda", "modulo a", 0);
clarg::argInt  modb("-modb", "modulo b", 0);
clarg::argInt  modc("-modc", "modulo c", 0);
clarg::argInt  nt("-nt", "number of threads", 0);
clarg::argInt  verb_level("-v", "verbosity level", 0);
clarg::argBool h("-h", "help message", false);

/* Run statistics. */
RunStatsTable compute_rst ("-comp_rdt", "Array multiply statistics raw data table");
RunStatsTable compute_rev_rst ("-comp_rev_rdt", "Array reverse multiply statistics raw data table");

struct thread_arg_t
{
  double* ara;
  double* arb;
  double* arc;
  unsigned sz;
};

void* mul_arr_thread(void* a)
{
  thread_arg_t* args = (thread_arg_t*) a;
  unsigned sz = args->sz;

  for (unsigned i=0; i<sz; i++)
    args->ara[i] = args->arb[i] * args->arc[i];
}

void mul_arr(double* ara, double* arb, double* arc, unsigned sz, unsigned nthreads)
{
  if (nthreads==0)
    nthreads = 1;

  pthread_t    *allthreads  = new pthread_t[nthreads];
  thread_arg_t *thread_args = new thread_arg_t[nthreads];

  unsigned chunk_sz = sz/nthreads;
  unsigned start = 0;

#define MIN(a,b) (a)<(b)?(a):(b)

  for (unsigned i=0; i<nthreads; i++, start+=chunk_sz) {
    thread_args[i].ara = &(ara[start]);
    thread_args[i].arb = &(arb[start]);
    thread_args[i].arc = &(arc[start]);
    thread_args[i].sz  = MIN(chunk_sz,sz-start);
  }

  /* Spaw threads */
  for(unsigned i=0; i<nthreads; i++) {
    pthread_create(&allthreads[i], NULL, mul_arr_thread, thread_args);
  }

  /* Join threads */
  for(unsigned i=0; i<nthreads; i++) {
    pthread_join(allthreads[i], NULL);
  }
}

void mul_arr_rev(double* ara, double* arb, double* arc, unsigned sz, unsigned threads)
{
  unsigned j=sz-1;
  for (unsigned i=0; i<sz; i++,j--)
    ara[i] = arb[i] * arc[j];
}



int main(int argc, char *argv[])
{
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
    
    unsigned sz = asz.get_value();

    double* ara = new double[sz];
    double* arb = new double[sz];
    double* arc = new double[sz];

    ElapsedTimeRunStat et_mul, et_rev_mul;

    mul_arr_rev(ara,arb,arc,sz,nt.get_value()); // Discard the first run
    compute_rev_rst.start();
    et_rev_mul.start();
    mul_arr_rev(ara,arb,arc,sz,nt.get_value());
    et_rev_mul.stop();
    compute_rev_rst.stop();

    mul_arr(ara,arb,arc,sz,nt.get_value()); // Discard the first run
    compute_rst.start();
    et_mul.start();
    mul_arr(ara,arb,arc,sz,nt.get_value());
    et_mul.stop();
    compute_rst.stop();

    #define GIGA 1000000000

    double mul_secs       = (et_mul.getElapsedMS()/1000.0);
    double mul_gflops     = (double) sz / (mul_secs * GIGA);
    double rev_mul_secs   = (et_rev_mul.getElapsedMS()/1000.0);
    double rev_mul_gflops = (double) sz / (rev_mul_secs * GIGA);
    cout << "a[1..." << sz << "] = b[1..." << sz << "] x c[1..." << sz << "] : " 
         << mul_secs << " : s : " << mul_gflops << " GFlops/s"<< endl;
    cout << "a[1..." << sz << "] = b[1..." << sz << "] x c[" << sz << "...1] : "
         << rev_mul_secs << " : s : " << rev_mul_gflops << " GFlops/s"<< endl;
    
    return 0; // Return ok
}

