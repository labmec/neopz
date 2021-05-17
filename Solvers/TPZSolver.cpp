#include "TPZSolver.h"
#include "Hash/TPZHash.h"

int TPZSolver::ClassId() const{
    return Hash("TPZSolver");
}