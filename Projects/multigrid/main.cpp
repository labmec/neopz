#include "pzmganalysis.h"
#include "pzcompel.h"

int main(){
  TPZCompEl::gOrder = 1;
   TPZMGAnalysis::main();
  return 0;
}
