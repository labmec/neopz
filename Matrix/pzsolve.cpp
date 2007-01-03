#include "pzsolve.h"

#include <stdlib.h>
using namespace std;

  /**
     Destructor
  */
TPZSolver::~TPZSolver() {}

TPZMatrixSolver::TPZMatrixSolver(TPZAutoPointer<TPZMatrix> Refmat) : fScratch() {
  fContainer = Refmat;
}

TPZMatrixSolver::TPZMatrixSolver() : fScratch() {
}

//misael
TPZMatrixSolver::TPZMatrixSolver(const TPZMatrixSolver &Source) : fScratch() {
    fReferenceMatrix = Source.fReferenceMatrix;
    fContainer = Source.fContainer;
}

// philippe 6/2/97
TPZMatrixSolver::~TPZMatrixSolver() {
}


void TPZMatrixSolver::SetMatrix(TPZAutoPointer<TPZMatrix> Refmat){
  fContainer = Refmat;  
}

void TPZMatrixSolver::ResetMatrix(){
  TPZAutoPointer<TPZMatrix> reset;
  fContainer = reset;
}


void TPZMatrixSolver::ShareMatrix(TPZMatrixSolver &other) {
    if(this == &other) return;
    fContainer = other.fContainer;
}
