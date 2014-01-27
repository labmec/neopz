/**
 * @author Gilvan S. Vieira
 */

#include <iostream> /* Standard IO */
#include "arglib.h" /* Command Line Arguments */
#include "run_stats_table.h"
#include <vector>

using namespace std;

#include "pzfmatrix.h" /* PZ Matrix Class */

#ifdef USING_TBB
#include "tbb/task_scheduler_init.h"
#include <tbb/tbb.h>
#include <tbb/parallel_for.h>
#include "tbb/partitioner.h"

using namespace tbb;
#endif

/**
 * @brief Command Line Arguments
 */
clarg::argBool h("-h", "Help Message", false);
clarg::argInt  dim("-d", "Dimension of the Matrices", 100);
clarg::argInt  num("-n", "Number of Matrices", 64);
clarg::argInt  repet("-r", "Times to Repeat Multiplication", 100);
clarg::argBool utbb("-tbb", "Using TBB", false);
clarg::argBool pcopy("-pcopy", "Using TBB Parallel Initialization", false);

void help(const char* prg)
{
    
} /* help */

TPZFMatrix<double> *randomMatrix()
{
    TPZFMatrix<double> *mMatrix = new TPZFMatrix<double>(dim.get_value(), dim.get_value());
    srand(time(NULL));
    for (int i=0; i<dim.get_value(); i++) {
        for (int j=0; j<dim.get_value(); j++) {
            (*mMatrix)(i,j) = rand() % 100;
        }
    }
    return mMatrix;
} /* randomMatrix */

void serialCopy(TPZFMatrix<double> *orig, vector<TPZFMatrix<double>* > &store)
{
    for (int w=0; w<num.get_value(); w++) {
        (*store[w]) = (*orig);
    }
} /* serialCopy */

void serialMultiplication(vector<TPZFMatrix<double>* > &input, vector<TPZFMatrix<double>* > &store) {
    for (int w=0; w<num.get_value(); w++) {
        (*store[w]) = (*input[w]) * (*input[w]);
    }
} /* serialMultiplication */

class Atribu {
public:
    TPZFMatrix<double> *input;
    vector<TPZFMatrix<double>* > *output;
#ifdef USING_TBB
    void operator()( const blocked_range<size_t>& r ) const {
        for (long i=r.begin(); i!=r.end(); ++i ) {
            (*(*output)[i]) = (*input);
        }
    }
#endif
};

class Multi {
public:
    vector<TPZFMatrix<double>* > *input;
    vector<TPZFMatrix<double>* > *output;
    
#ifdef USING_TBB
    void operator()( const blocked_range<size_t>& r ) const {
        for (long i=r.begin(); i!=r.end(); ++i ) {
            (*(*output)[i]) = (*(*input)[i]) * (*(*input)[i]);
        }
    }
#endif
};

RunStatsTable mult_rdt   ("-mult_rdt", "Multiplication statistics raw data table");

int main(int argc, char **argv)
{
    /* parse the arguments */
    if (clarg::parse_arguments(argc, argv)) {
        cerr << "Error when parsing the arguments!" << endl;
        return 1;
    }
    /* checking if the Help Message was request */
    if (h.get_value() == true) {
        help(argv[0]);
        return 1;
    }
    
    TPZFMatrix<double> *m = randomMatrix();
    
    vector<TPZFMatrix<double>* > matrices(num.get_value(), new TPZFMatrix<double>(dim.get_value(), dim.get_value()));
    vector<TPZFMatrix<double>* > results(num.get_value(), new TPZFMatrix<double>(dim.get_value(), dim.get_value()));
    
#ifdef USING_TBB
    affinity_partitioner ap;
#endif
    if (pcopy.get_value() == false) {
        serialCopy(m, matrices);
    } else {
        Atribu builder;
        builder.input = m;
        builder.output = &matrices;
        
#ifdef USING_TBB
        parallel_for(blocked_range<size_t>(0,num.get_value(), 1), builder, ap);
#endif
    }
    
    mult_rdt.start();
    if(utbb.get_value() == false) {
        
        for (int i=0; i<repet.get_value();i++)
            serialMultiplication(matrices, results);
        
    } else {
#ifdef USING_TBB
        task_scheduler_init init;
        
        Multi op;
        op.input = &matrices;
        op.output = &results;
        
        for (int i=0; i<repet.get_value();i++)
            parallel_for(blocked_range<size_t>(0,num.get_value(), 1), op, ap);
#endif
    }
    
    mult_rdt.stop();
    
} /* main */