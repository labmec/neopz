#include "pzsolve.h"

#include <stdlib.h>

  /**
     Destructor
  */
TPZSolver::~TPZSolver() {}

TPZMatrixSolver::TPZMatrixSolver(TPZMatrix *Refmat) : fScratch() {

  fContainer = new TPZContainer(Refmat);
  gnumcreated++;
}
//misael
TPZMatrixSolver::TPZMatrixSolver(const TPZMatrixSolver &Source) : fScratch() {
    fContainer = Source.fContainer;
    fContainer->IncreaseRefCount();
    gnumcreated++;

}

// philippe 6/2/97
TPZMatrixSolver::~TPZMatrixSolver() {
    fContainer->DecreaseRefCount();
    gnumdeleted++;
}

int TPZMatrixSolver::gnumcreated = 0;
int TPZMatrixSolver::gnumdeleted = 0;
//misael

void TPZMatrixSolver::SetMatrix(TPZMatrix *Refmat){
    if(fContainer->Matrix() == Refmat || !fContainer->Matrix()) {
        fContainer->SetMatrix(Refmat);
    } else {
        fContainer->DecreaseRefCount();
        fContainer = new TPZContainer(Refmat);
    }
}

void TPZMatrixSolver::ResetMatrix(){
  fContainer->SetMatrix(0);
}

TPZMatrixSolver::TPZContainer::TPZContainer(TPZMatrix *mat){
    fRefMat = mat;
    fRefCount = 1;
    gnumcreated++;
}

TPZMatrixSolver::TPZContainer::~TPZContainer(){
    if(fRefMat) delete fRefMat;
    gnumdeleted++;
}

int TPZMatrixSolver::TPZContainer::gnumcreated = 0;
int TPZMatrixSolver::TPZContainer::gnumdeleted = 0;

void TPZMatrixSolver::TPZContainer::Diagnose(ostream &out) {
  out << "TPZMatrixSolver::TPZContainer\nnumber of objects created " << gnumcreated <<
    "\nnumber of objects deleted " << gnumdeleted << endl;
}

void TPZMatrixSolver::TPZContainer::IncreaseRefCount(){
    fRefCount++;
}

void TPZMatrixSolver::TPZContainer::DecreaseRefCount(){
    fRefCount--;
    if(fRefCount < 0) {
        cout << "TPZMatrixSolver::TPZContainer::DecreaseRefCount wrong data structure\n";
        return;
    }
    if(fRefCount == 0) delete this;
}

TPZMatrix *TPZMatrixSolver::TPZContainer::Matrix() {
    return fRefMat;
}

void TPZMatrixSolver::TPZContainer::SetMatrix(TPZMatrix *mat) {
    fRefMat = mat;
}

void TPZMatrixSolver::ShareMatrix(TPZMatrixSolver &other) {
    if(this == &other) return;
    fContainer->DecreaseRefCount();
    fContainer = other.fContainer;
    fContainer->IncreaseRefCount();
}
