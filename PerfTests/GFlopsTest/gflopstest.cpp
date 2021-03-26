/**
 * @file
 * @brief Measure the amount of GFlops/s the system can sustain.
 * @author Edson Borin
 * @since 2012
 */

#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <thread>
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
clarg::argInt  cm("-cm", "clean memory before execution", 512);


void clean_mem(unsigned sz)
{
	unsigned* buffer = (unsigned*) malloc(sz/4);
	for (unsigned i=0; i<sz/4; i++)
		buffer[i] = i;
	free(buffer);
} 

/* Run statistics. */
RunStatsTable mul_rst ("-mul_rdt", "Array multiply statistics raw data table");
RunStatsTable imul_rst ("-imul_rdt", "Array immeditate multiply statistics raw data table");
RunStatsTable mulred_rst ("-mulred_rdt", "Array multiply and reduce statistics raw data table");
RunStatsTable add_rst ("-add_rdt", "Array add statistics raw data table");
RunStatsTable acc_rst ("-acc_rdt", "Array accumulate statistics raw data table");
RunStatsTable mulsingle_rst ("-mulsingle_rdt", "Fake array-multiply statistics raw data table");

struct thread_arg_t
{
	double* ara;
	double* arb;
	double* arc;
	unsigned sz;
};

void process_arr(double* ara, double* arb, double* arc, unsigned sz, unsigned nthreads, void* (*fun)(void*))
{
	if (nthreads==0)
		nthreads = 1;
	
	std::vector<std::thread> allthreads;
	thread_arg_t *thread_args = new thread_arg_t[nthreads];
	
	unsigned chunk_sz = sz/nthreads;
	unsigned start = 0;
	
#define MIN_T(a,b) (a)<(b)?(a):(b)
	
	for (unsigned i=0; i<nthreads; i++, start+=chunk_sz) {
		thread_args[i].ara = &(ara[start]);
		thread_args[i].arb = &(arb[start]);
		thread_args[i].arc = &(arc[start]);
		thread_args[i].sz  = MIN_T(chunk_sz,sz-start);
	}
	
	/* Spaw threads */
	for(unsigned i=0; i<nthreads; i++) {
    allthreads.push_back(std::thread(fun,&(thread_args[i])));
	}
	
	/* Join threads */
	for(unsigned i=0; i<nthreads; i++) {
    allthreads[i].join();
	}
}

template<class T>
struct thread_map1_arg
{
	T* array;
	unsigned sub_sz;
	T (*map_func)(T); 
};

template<class T>
void* thread_map1_worker(void* arg)
{
	thread_map1_arg<T>* args = (thread_map1_arg<T>*) arg;
	T* array = args->array;
	unsigned sub_sz = args->sub_sz;
	
	for (unsigned i=0; i<sub_sz; i++)
		array[i] = args->map_func(array[i]);
	return 0;
}

template<class T>
void thread_map1(T* array, T (*map1)(T), unsigned sz, unsigned nthreads)
{
	if (nthreads==0)
		nthreads = 1;
	
	std::vector<std::thread> allthreads;
	thread_map1_arg<T> *thread_args = new thread_map1_arg<T>[nthreads];
	
	unsigned chunk_sz = sz/nthreads;
	unsigned start = 0;
	
	for (unsigned i=0; i<nthreads; i++, start+=chunk_sz) {
		thread_args[i].array = &(array[start]);
		thread_args[i].sub_sz = MIN_T(chunk_sz,sz-start);
		thread_args[i].map_func = map1;
	}
	
	/* Spaw threads */
	for(unsigned i=0; i<nthreads; i++) {
    allthreads.push_back(std::thread(thread_map1_worker<T>,&(thread_args[i])));
	}
	
	/* Join threads */
	for(unsigned i=0; i<nthreads; i++) {
    allthreads[i].join();
	}
}


void mul_arr_rev(double* ara, double* arb, double* arc, unsigned sz, unsigned threads)
{
	unsigned j=sz-1;
	for (unsigned i=0; i<sz; i++,j--)
		ara[i] = arb[i] * arc[j];
}

double sqrt_dbl(double v) { return v*v; }

void* mul_arr(void* a)
{
	thread_arg_t* args = (thread_arg_t*) a;
	unsigned sz = args->sz;
	double* ara = args->ara;
	double* arb = args->arb;
	double* arc = args->arc;
	for (unsigned i=0; i<sz; i++)
		ara[i] = arb[i] * arc[i];
	
	return 0;
}

void* imul_arr(void* a)
{
	thread_arg_t* args = (thread_arg_t*) a;
	unsigned sz = args->sz;
	double* ara = args->ara;
	
	for (unsigned i=0; i<sz; i++)
		ara[i] = ara[i] * 1.9752; // Multiply by a constant
	
	return 0;
}

void* add_arr(void* a)
{
	thread_arg_t* args = (thread_arg_t*) a;
	unsigned sz = args->sz;
	double* ara = args->ara;
	double* arb = args->arb;
	double* arc = args->arc;
	for (unsigned i=0; i<sz; i++)
		ara[i] = arb[i] + arc[i];
	
	return 0;
}

void* acc_arr(void* a)
{
	thread_arg_t* args = (thread_arg_t*) a;
	unsigned sz = args->sz;
	double* ara = args->ara;
	double* arb = args->arb;
	for (unsigned i=0; i<sz; i++)
		ara[i] = ara[i] + arb[i];
	
	return 0;
}

double global_res;
void* mulred_arr(void* a)
{
	thread_arg_t* args = (thread_arg_t*) a;
	unsigned sz = args->sz;
	double* arb = args->arb;
	double* arc = args->arc;
	double res = 0;
	for (unsigned i=0; i<sz; i++)
		res += arb[i] * arc[i];
	
	global_res = res;
	return 0;
}

void* mulsingle_arr(void* a)
{
	thread_arg_t* args = (thread_arg_t*) a;
	unsigned sz = args->sz;
	double* arb = args->arb;
	double* arc = args->arc;
	double res = 0;
	for (unsigned i=0; i<sz; i++)
		res += arb[0] * arc[0];
	
	global_res = res;
	return 0;
}

void profile(double* ara, double* arb, double* arc, unsigned sz, unsigned num_threads, 
             void* (*fun)(void*), ElapsedTimeRunStat& et, RunStatsTable& rst)
{
	process_arr(ara,arb,arc,sz,num_threads,fun);
	rst.start();
	et.start();
	process_arr(ara,arb,arc,sz,num_threads,fun);
	et.stop();
	rst.stop();
}

double gflops(ElapsedTimeRunStat& et, unsigned sz)
{
#define GIGA 1000000000
	
	double secs   = (et.getElapsedMS()/1000.0);
	double gflops = (double) sz / (secs * GIGA);
	return gflops;
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
	
    if (cm.was_set()) {
#define MEGABYTE (1024*1024)
		unsigned sz = cm.get_value() * MEGABYTE;
		cout << "Cleaning memory: " << cm.get_value() << " MB ...";
		clean_mem(sz);
		cout << "[Done]" << endl;
    }
    
    unsigned sz = asz.get_value();
	
    /* Create arrays. */
    double* ara = new double[sz];
    double* arb = new double[sz];
    double* arc = new double[sz];
	
    ElapsedTimeRunStat et_mul, et_imul, et_add, et_acc, et_mulred, et_mulsingle;
	
    profile(ara,arb,arc,sz,nt.get_value(),
            mul_arr,et_mul,mul_rst);
	
    profile(ara,arb,arc,sz,nt.get_value(),
            imul_arr,et_imul,imul_rst);
	
    profile(ara,arb,arc,sz,nt.get_value(),
            add_arr,et_add,add_rst);
	
    profile(ara,arb,arc,sz,nt.get_value(),
            acc_arr,et_acc,acc_rst);
	
    profile(ara,arb,arc,sz,nt.get_value(),
            mulred_arr,et_mulred,mulred_rst);
	
    profile(ara,arb,arc,sz,nt.get_value(),
            mulsingle_arr,et_mulsingle,mulsingle_rst);
	
    //compute_rev_rst.start();
    //et_rev_mul.start();
    //thread_map1(ara, sqrt_dbl, sz, nt.get_value());
    //et_rev_mul.stop();
    //compute_rev_rst.stop();
	
	
    cout << "Array mul performance     : " << gflops(et_mul,sz) << endl;
    cout << "Array imul performance    : " << gflops(et_imul,sz) << endl;
    cout << "Array add performance     : " << gflops(et_add,sz) << endl;
    cout << "Array acc performance     : " << gflops(et_acc,sz) << endl;
    cout << "Array mul red performance : " << gflops(et_mulred,sz) << endl;
    cout << "Fake array mul performance: " << gflops(et_mulsingle,sz) << endl;
	
    return 0; // Return ok
}

