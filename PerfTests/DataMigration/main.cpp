#include <iostream>
#include "arglib.h"
#include "pzfmatrix.h"
#include "pzparallel.h"

clarg::argInt number_of_matrices("-n", "Number of matrices.", 1);
clarg::argInt dimension("-d", "Matrices dimension M x M", 2500);

class Benchmark {
private:
    std::vector<TPZFMatrix<REAL> *> matrices;
public:
    Benchmark(int number_of_matrices, int dimension) : matrices(number_of_matrices) {
        for (int i = 0; i < number_of_matrices; i++) {
            matrices[i] = new TPZFMatrix<REAL>(dimension, dimension, M_PI);
        }
    }
    
    void operator()(int i) {
        matrices[i]->Decompose_LU();
    }
    
    ~Benchmark() {
        for (int i = 0; i < matrices.size(); i++)
            delete matrices[i];
        matrices.clear();
    }
};
int main(int argc, char **argv) {
    // parse command line arguments
    clarg::parse_arguments(argc, argv);
    // define input parameters

    Benchmark bt(number_of_matrices.get_value(), dimension.get_value());
    
    pz::parallel_for(number_of_matrices.get_value(), bt);

}