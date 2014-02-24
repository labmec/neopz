#include <iostream>         // input and output
#include <vector>           // standard vector container

#include "pzskylmat.h"      // skyline matrix
#include "input.h"          // CreateCuboSkyMatrix
#include "arglib.h"         // arguments lib

#ifdef USING_TBB
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/task_scheduler_init.h>
#endif

using namespace std;

clarg::argString    filename("-f", "Cube input file name.", "cube1.txt");
clarg::argInt       plevel("-p", "Level of order.", 2);
clarg::argInt       nmatrices("-n", "Number of matrices.", 64);
clarg::argInt       nloop("-l", "Number of loop iterations of the Subst_Backward/Subst_Forward", 1);
clarg::argBool      usetbb("-tbb", "Use of parallel tbb version", false);

#ifdef USING_TBB
class tbb_cholesky {
public:
    
    vector<TPZSkylMatrix<REAL>* > *fTasks;
    
    void operator()(const tbb::blocked_range<size_t>& range) const {
        for(size_t i=range.begin(); i!=range.end(); ++i ) {
            (*fTasks)[i]->Decompose_Cholesky();
        }
    } // operator()
};

class tbb_substitution {
public:
    
    vector<TPZSkylMatrix<REAL>* > *fTasks;
    
    void operator()(const tbb::blocked_range<size_t>& range) const {
        for(size_t i=range.begin(); i!=range.end(); ++i ) {
            TPZFMatrix<REAL> f((*fTasks)[i]->Dim(),1,M_PI);
            
            (*fTasks)[i]->Subst_Backward(&f);
            (*fTasks)[i]->Subst_Forward(&f);
        }
    } // operator()
};
#endif

int main(int argc, char **argv)
{
    // parse the arguments
    if (clarg::parse_arguments(argc, argv)) {
        cerr << "Error when parsing the arguments!" << endl;
        return 1;
    }
    
#ifdef USING_TBB
    tbb::task_scheduler_init init;
#endif
    
    // generate original SkyMatrix
    TPZAutoPointer<TPZMatrix<REAL> > orig = Input::CreateCuboSkyMatrix(filename.get_value(), plevel.get_value());
    TPZSkylMatrix<REAL> *sky = dynamic_cast<TPZSkylMatrix<REAL> *> (orig.operator->());
    
    vector<TPZSkylMatrix<REAL>* > fTasks(nmatrices.get_value());
    
    // serial copy
    for (int i=0; i<nmatrices.get_value(); i++) {
        fTasks[i] = new TPZSkylMatrix<REAL>(*sky);
    }
    
    if (!usetbb.get_value()) {
        // serial decompose cholesky
        cout << "----> Decompose_Cholesky" << endl;
        for (int i=0; i<nmatrices.get_value(); i++) {
            fTasks[i]->Decompose_Cholesky();
        }
        // serial Subst_Backward/Subst_Forward
        cout << "----> Subst_Backward/Subst_Forward" << endl;
        for (int k=0; k<nloop.get_value();k++) {
            for (int i=0; i<nmatrices.get_value(); i++) {
                TPZFMatrix<REAL> f(fTasks[i]->Dim(),1,M_PI);
                fTasks[i]->Subst_Backward(&f);
                fTasks[i]->Subst_Forward(&f);
            }
        }
    } else {
#ifdef USING_TBB
        
        tbb_cholesky disp;
        disp.fTasks = &fTasks;
        
        tbb::affinity_partitioner ap;
        cout << "----> Decompose_Cholesky" << endl;
        parallel_for(tbb::blocked_range<size_t>(0, nmatrices.get_value()), disp, ap);
        
        tbb_substitution dispb;
        dispb.fTasks = &fTasks;
        cout << "----> Subst_Backward/Subst_Forward" << endl;
        for (int k=0; k<nloop.get_value();k++) {
            parallel_for(tbb::blocked_range<size_t>(0, nmatrices.get_value()), dispb, ap);
        }
#else
        cout << "Compiled without TBB support." << endl;
#endif
    }
    
} // main
