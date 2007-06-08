
#include "pzbdstrmatrix.h"
#include "pzmganalysis.h"
#include "pzcompel.h"
#include "pzonedref.h"

int main() {

//    TPZCompEl::gOrder = 3;
  cmesh.SetDefaultOrder(3);
  //  return TPZBlockDiagonalStructMatrix::main();
//     return TPZOneDRef::main();
    return TPZMGAnalysis::main();

}
