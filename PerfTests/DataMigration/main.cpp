#include <iostream>
#include "arglib.h"
#include "pzfmatrix.h"
#include "pzparallel.h"
#include "run_stats_table.h"

#include <float.h> // FLT_MAX
#include <stdio.h> // printf

clarg::argInt number_of_matrices("-n", "Number of matrices.", 2);
clarg::argInt dimension("-d", "Matrices dimension M x M", 1000);

#ifndef DATATYPE
#define RPT 5
#endif

RunStatsTable total_rst   ("-tot_rdt", "Whole program (total) statistics raw data table");

class MatrixBenchmark {
private:
    std::vector<TPZFMatrix<REAL> *> matrices;
public:
    MatrixBenchmark(int number_of_matrices, int dimension) : matrices(number_of_matrices) {
        for (int i = 0; i < number_of_matrices; i++) {
            matrices[i] = new TPZFMatrix<REAL>(dimension, dimension, M_PI);
        }
    }
    
    void operator()(int i) {
        matrices[i]->Decompose_LU();
    }
    
    ~MatrixBenchmark() {
        for (int i = 0; i < matrices.size(); i++)
            delete matrices[i];
        matrices.clear();
    }
};

/* Code to read the wall clock time.              */
#include <sys/time.h>
double mysecond()
{
    struct timeval tp;
    struct timezone tzp;
    gettimeofday(&tp,&tzp);
    return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}

/* Kernel name. */
const char* kernel_name = "TPZFMatrix LU Decomposition";

void kernel()
{
    MatrixBenchmark bt(number_of_matrices.get_value(), dimension.get_value());
    pz::parallel_for(number_of_matrices.get_value(), bt);
}

/* -----------------------------*/
int main()
{
    uint64_t k;
    double             times[RPT];
    double             mintime = FLT_MAX;
    double             avgtime = 0;
    double             maxtime = 0;
    double             t;
    
    printf("Kernel name     : %s\n",kernel_name);
    printf("# of runs       : %d\n", RPT);
    
    /* Main loop. */
    for (k=0; k<RPT; k++)
    {
        t = mysecond();
        /* Kernel */
        kernel();
        times[k] = mysecond() - t;
        //printf(" -> %6.2f s\n", times[k]);
    }
    
    /* Final report */
    for (k=0; k<RPT; k++)
    /* Discard first iteration (k=1). */
    {
        avgtime = avgtime + times[k];
        mintime = MIN(mintime, times[k]);
        maxtime = MAX(maxtime, times[k]);
    }
    /* Print Report */
    avgtime = avgtime / (RPT-1);
    printf("Avg time        : %6.2f\n",avgtime);
    printf("Min time        : %6.2f\n",mintime);
    printf("Max time        : %6.2f\n",maxtime);
    
    return 0;
}

